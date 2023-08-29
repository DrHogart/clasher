from datetime import datetime
from Bio.Seq import Seq
from Bio import SeqIO, SeqRecord
from Bio.SeqRecord import SeqRecord
from pybedtools import BedTool, example_filename, cbedtools
import copy
import pandas as pd
import os
import re
import subprocess
import statistics
import sys

def index_exists (index):
    index_folder = '/'.join(index.split('/')[:-1])
    index_base = index.split('/')[-1]
    if not os.path.isdir(index_folder):
        return False
    
    files = os.listdir(index_folder)
    for file in files:
        if file.startswith(index_base):
            return True
    return False

def bedtool_from_intervals (intervals):
    fields_count = min(len(intervals[0].fields), 6)
    columns = ['chrom', 'start', 'end', 'name', 'score', 'strand'][:fields_count]
    df = pd.DataFrame([i.fields[:fields_count] for i in intervals], columns=columns)
    return BedTool.from_dataframe(df)

def expand_intervals (bed, scope, expand_to_scope=False):
    intervals_expanded = []
    for c in bed:
        minus = 0
        plus = 0
        if expand_to_scope:
            scope_len_diff = max(0, (scope - c.length))
            minus = scope_len_diff // 2
            plus = scope_len_diff // 2 + scope_len_diff % 2
        else:
            minus = scope // 2
            plus = scope - minus
        fields = c.fields[:6]
        fields[1], fields[2] = int(fields[1]) - minus, int(fields[2]) + plus
        intervals_expanded.append(cbedtools.Interval(*fields))
    return bedtool_from_intervals(intervals_expanded)

def parse_coservation (file_name):
    file = open(file_name, 'r')
    lines = file.readlines()
    file.close()
    
    # Get only the fourth column
    lines = [line.split('\t')[3].strip() for line in lines]
    cons_dict = dict()
    
    # Iterate through file lines
    for line in lines:
        peak = line.split('-')[0]
        conservation = line[line.rfind('|') + 1:]
        conservation = float(conservation) \
            if conservation.replace('.', '').replace('-', '').isnumeric() else ''
        
        # Add conservation value to peak in dict
        if peak in cons_dict.keys():
            cons_dict[peak].append(conservation)
        else:
            cons_dict[peak] = [conservation]
    return cons_dict

# INPUT VARIABLES START
# =====================

research_id = 'GSE73059_cw_pureclip'
fastq_folder = '/home/dimasnitkin/workspace/human/ncbi_geo/GSE73059_cw'
acc_numbers = 'SRR2413179,SRR2413181,SRR2413180,SRR2413182,SRR2413281,SRR2413282,SRR2413284,SRR2413285,SRR2413286,SRR2413287'
adapter_5_prime = 'AGGGAGGACGATGCGG'
adapter_3_prime = 'GTGTCAGTCACTTCCAGCGG'
bowtie_folder = '/home/dimasnitkin/workspace/human/bowtie'
index_folder = '{}/refs'.format(bowtie_folder)
rrna_index = '{}/refs/rrna_rnacentral/rrna'.format(bowtie_folder)
bam_folder = '{}/results'.format(bowtie_folder)

# Gencode Genome
genome_index = '{}/refs/gencode_genome/'.format(bowtie_folder)
genome_index += 'GRCh38.p13.genome.chr/GRCh38.p13.genome.chr'
genome_folder = '/home/dimasnitkin/workspace/human/gencode'
genome_ann = '{}/gencode.v43.annotation.gff3'.format(genome_folder)
genome_fasta = '{}/GRCh38.primary_assembly.genome.fa'.format(genome_folder)

log_folder = '/home/dimasnitkin/workspace/logs'
log_file = '{}/{}.log'.format(log_folder, research_id)
mirna_file = '/home/dimasnitkin/workspace/human/mirbase/mature_mirnas_bt.fa'
mirna_hairpins_index = '/home/dimasnitkin/workspace/human/bowtie/refs/mirnas_mirbase_hairpins_bt/mirnas'
csv_folder = '/home/dimasnitkin/workspace/human/csv'
csv_file_base = '{}/{}'.format(csv_folder, research_id)

# Conda environments
clipper_env = 'clipper3'
clipper_species = 'gencode.v43'

# Folder with wigFix chromosomes files converted to BED 
phylop_folder = '/home/dimasnitkin/workspace/human/phast'

# INPUT VARIABLES END
# ===================

if not os.path.isdir(bam_folder):
    os.mkdir(bam_folder)
if not os.path.isdir(csv_folder):
    os.mkdir(csv_folder)
if not os.path.isdir(index_folder):
    os.mkdir(index_folder)
if not os.path.isdir(log_folder):
    os.mkdir(log_folder)
if not os.path.isdir(log_folder + '/fastqc'):
    os.mkdir(log_folder + '/fastqc')
if not os.path.isfile(genome_fasta + '.fai'):
    os.system('samtools faidx {}'.format(genome_fasta))
if not os.path.isfile(genome_fasta + '.sizes'):
    os.system('cut -f 1,2 {}.fai > {}.sizes'.format(genome_fasta, genome_fasta))

# Asserts
assert os.path.isdir(fastq_folder), 'FASTQ folder not found: {}'.format(fastq_folder)
assert len(os.listdir(fastq_folder)) != 0, 'FASTQ folder is empty: {}'.format(fastq_folder)
# for acc in acc_numbers.split(','):
#     fastq_file = '{}/{}.fastq'.format(fastq_folder, acc)
#     assert os.path.isfile(fastq_file), 'FASTQ file not found: {}'.format(fastq_file)
assert os.path.isfile(mirna_file), 'miRNAs file not found: {}'.format(mirna_file)
assert index_exists(rrna_index), 'rRNA index not found: {}'.format(rrna_index)
assert index_exists(genome_index), 'Genome index not found: {}'.format(genome_index)
assert os.path.isfile(genome_ann), 'Genome annotation not found: {}'.format(genome_ann)
conda_envs = subprocess.check_output(['conda', 'env', 'list'], encoding='UTF-8')
assert '\n' + clipper_env + ' ' in conda_envs, 'Conda environment `{}` not found.'.format(clipper_env)

# Reads preprocessing
os.system('/bin/bash preprocessing.sh -f {} -r {} -o {} -i {} -b {} -c {} -a {} -m {} -w {} -l {}'.format(
    adapter_5_prime, adapter_3_prime, fastq_folder, index_folder, bam_folder, rrna_index, acc_numbers, mirna_file, research_id, log_file))

# Removing duplicates
print('Removing duplicates...')
for acc in acc_numbers.split(','):
    print(acc)
    bed_file = '{}/results/{}.bed'.format(bowtie_folder, acc)
    names = ['chimera_id', 'start', 'end', 'mirna', 'score', 'strand']

    bed_df = pd.read_csv(bed_file, header=None, names=names, sep='\t')
    print('Before dropping duplicates:', len(bed_df))
    bed_df.drop_duplicates(subset=['chimera_id', 'start', 'end', 'strand'], keep='first', inplace=True)
    bed_df.drop_duplicates(subset=['chimera_id'], keep=False, inplace=True)
    bed_df.to_csv(bed_file, sep='\t', index=False, header=False)
    print('After dropping duplicates:', len(bed_df))

    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    os.system('echo "(((((Chimeras with 1 miRNA|||||{})))))\n{}" >> {}'.format(
        now, bed_df.shape[0], log_file))
print('Removing duplicates complete!')

# Creating pandas dataframes
print('Creating pandas dataframes...')
mirnas = SeqIO.to_dict(SeqIO.parse(mirna_file, 'fasta'))
aln_counter = 0

for acc in acc_numbers.split(','):
    print(acc)
    chimeras_file = '{}/{}.norrna.fastq'.format(fastq_folder, acc)
    chimeras = SeqIO.to_dict(SeqIO.parse(chimeras_file, 'fastq'))
    for key, value in chimeras.items():
        value.description = ''
    
    bed_file = '{}/{}.bed'.format(bam_folder, acc)
    bed = BedTool(bed_file)
    df = pd.DataFrame(columns=['alignment_id', 'chimera_id', 'is_revcomp',
                              'mirna_id', 'mirna_start', 'mirna_end'])
    
    no_targets_counter = 0
    target_sites = []
    
    for c in bed:
        chimera_id = c.fields[0]
        
        if chimera_id not in chimeras.keys():
            continue
            
        strand = c.fields[5]
        is_revcomp = True if strand == '-' else False
        part_length = 15

        chimera, mirna_start, mirna_end = None, None, None

        mirna_name = c.fields[3]
        mirna = mirnas[mirna_name]

        if is_revcomp:
            chimera = copy.deepcopy(chimeras[chimera_id].reverse_complement())
            chimera.id = chimera_id
            chimera.description = 'revcomp'
            mirna_start = len(chimera) - len(mirna) - int(c.fields[1])
        else:
            mirna_start = int(c.fields[1])
            chimera = copy.deepcopy(chimeras[chimera_id])

        mirna_end = mirna_start + len(mirna) - 1

        chimera_seq = str(chimera.seq)

        len_diff = len(chimera) - len(mirna)
        if len_diff < part_length:
            continue

        mirna_part = str(chimera[mirna_start:mirna_end + 1].seq)

        left_part = copy.deepcopy(chimera[0:mirna_start])
        right_part = copy.deepcopy(chimera[mirna_end + 1:])
        
        if len(left_part) < part_length and len(right_part) < part_length:
            no_targets_counter += 1
            continue
        
        if len(left_part) >= part_length:
            left_part.id = f"{left_part.id}_aln{aln_counter}_left"
            target_sites.append(left_part)

        if len(right_part) >= part_length:
            right_part.id = f"{right_part.id}_aln{aln_counter}_right"
            target_sites.append(right_part)

        if len(left_part) >= part_length or len(right_part) >= part_length:
            row = [f"aln{aln_counter}",
                  chimera_id,
                  is_revcomp,
                  mirna_name,
                  mirna_start,
                  mirna_end]
            df.loc[len(df)] = row
            
        aln_counter += 1
        del chimeras[chimera_id]
        
    target_sites_file = '{}/{}.T.fastq'.format(fastq_folder, acc)
    with open(target_sites_file, 'w') as output_handle:
        SeqIO.write(target_sites, output_handle, 'fastq')
    
    log_file = "{}/{}.log".format(log_folder, acc)
    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    os.system('echo "(((((Chimeras with too short non-miRNA sites|||||{})))))\n{}" >> {}'.format(
        now, no_targets_counter, log_file))
    
    csv_file = '{}/{}.1.csv'.format(csv_folder, acc)
    df.to_csv(csv_file, index=False, sep='\t')
    print(csv_file, 'file created.')
    
del chimeras
print('Creating pandas dataframes complete!')
    
# Annotating target sites (T)
os.system('/bin/bash target_sites_annotation.sh -a {} -i {} -b {} -l {} -g {} -f {} -m {}'.format(
    acc_numbers, fastq_folder, bam_folder, log_file, genome_index, genome_ann, mirna_hairpins_index))

# Adding features to dataframes
print('Adding features to dataframes...')

for acc in acc_numbers.split(','):
    print(acc)
    df = pd.read_csv('{}/{}.1.csv'.format(csv_folder, acc), sep='\t')
    bed = {}

    for feature in ['CDS', 'intron', 'three_prime_UTR', 'five_prime_UTR', 'nc_RNA', 'miRNA', 'other']:
        file = '{}/{}.T.{}.bed'.format(bam_folder, acc, feature)
        if os.path.exists(file):
            bed[feature] = BedTool(file)

    df.insert(len(df.columns), 'is_mirna_first', None)
    df.insert(len(df.columns), 'feature_left', None)
    df.insert(len(df.columns), 'feature_right', None)
    df.insert(len(df.columns), 'feature', None)
    df.insert(len(df.columns), 'target_id', None)
    df.insert(len(df.columns), 'gene_id', None)
    df.insert(len(df.columns), 'gene_type', None)
    df.insert(len(df.columns), 'chromosome', None)
    df.insert(len(df.columns), 'chrom_strand', None)
    df.insert(len(df.columns), 'chrom_start', None)
    df.insert(len(df.columns), 'chrom_end', None)

    counter = 0
    # Add features
    print('Adding features...')
    for bed_key, bed_value in bed.items():
        for c in bed_value:
            target_site_id = c.name

            side = target_site_id.split('_')[-1]
            is_mirna_first = True if side == 'right' else False
            alignment_id = re.search('aln[0-9]+', target_site_id).group(0)

            target_id = c.fields[20].split('Parent=')[1].split(';')[0] \
                if 'Parent=' in c.fields[20] else ''
            gene_id = c.fields[20].split('gene_id=')[1].split(';')[0]

            df.loc[df['alignment_id'] == alignment_id, 'is_mirna_first'] = is_mirna_first
            df.loc[df['alignment_id'] == alignment_id, f"feature_{side}"] = bed_key
            df.loc[df['alignment_id'] == alignment_id, 'feature'] = bed_key

            counter += 1

    # Leave alignments with one feature: left or right
    print('Leaving alignments with one feature: left or right...')
    df = df[df['feature_right'].notna() != df['feature_left'].notna()]
    df = df.drop(['feature_right', 'feature_left'], axis='columns')

    counter = 0
    # Add alignment characteristics
    print('Adding alignment characteristics...')
    for bed_key, bed_value in bed.items():
        for c in bed_value:

            target_site_id = c.name
            alignment_id = re.search('aln[0-9]+', target_site_id).group(0)
            target_id = c.fields[20].split('Parent=')[1].split(';')[0] \
                if 'Parent=' in c.fields[20] else None
            gene_id = c.fields[20].split('gene_id=')[1].split(';')[0]
            gene_type = c.fields[20].split('gene_type=')[1].split(';')[0] \
                if 'gene_type=' in c.fields[20] else None

            df.loc[df['alignment_id'] == alignment_id, 'target_id'] = target_id
            df.loc[df['alignment_id'] == alignment_id, 'gene_id'] = gene_id
            df.loc[df['alignment_id'] == alignment_id, 'gene_type'] = gene_type
            df.loc[df['alignment_id'] == alignment_id, 'chromosome'] = c.chrom
            df.loc[df['alignment_id'] == alignment_id, 'chrom_strand'] = c.strand
            df.loc[df['alignment_id'] == alignment_id, 'chrom_start'] = c.start
            df.loc[df['alignment_id'] == alignment_id, 'chrom_end'] = c.end

            counter += 1

    file_base = '{}/{}.T.filtered'.format(bam_folder, acc)
    bed_filtered_file = open(file_base + '.bed', 'w')
    for index, row in df.iterrows():
        bed_filtered_file.write('{}\t{}\t{}\t{}\t.\t{}\n'.format(
            row['chromosome'], row['chrom_start'], row['chrom_end'], row['alignment_id'], row['chrom_strand']))
    bed_filtered_file.close()

    now = datetime.now().strftime('%Y-%m-%d %H:%M')
    os.system('echo "(((((Chimeras after annotating target sites|||||{})))))\n{}" >> {}'.format(
        now, len(df), log_file))

    df.to_csv('{}/{}.2.csv'.format(csv_folder, acc), index=False, sep='\t')
print('Adding features to dataframes complete!')

# Converting BED files to BAM file
print('Converting BED files to BAM file...')
bed_files = ['{}/{}.T.filtered.bed'.format(bam_folder, acc) for acc in acc_numbers.split(',')]
bed_files = ' '.join(bed_files)
bam_base = '{}/{}.T.filtered'.format(bam_folder, research_id)
os.system('cat {} > {}.bed'.format(bed_files, bam_base))
os.system('bedtools bedtobam -i {}.bed -g {}.sizes > {}.bam'.format(
    bam_base, genome_fasta, bam_base))
bam_sorted = bam_base + '.sorted.bam'
os.system('samtools sort -m 2G {}.bam -o {}'.format(
    bam_base, bam_sorted))
os.system('samtools index ' + bam_sorted)
print('Converting BED files to BAM file complete!')

# Merging dataframes
print('Merging dataframes...')
first_acc = acc_numbers.split(',')[0]
os.system('cp {}/{}.2.csv {}/{}.2.csv'.format(csv_folder, first_acc, csv_folder, research_id))
for acc in acc_numbers.split(',')[1:]:
    os.system('tail -n +2 {}/{}.2.csv >> {}/{}.2.csv'.format(csv_folder, acc, csv_folder, research_id))
print('Merging dataframes complete!')

# Deleting intermediate files
print('Deleting intermediate files...')
features_deleted = ['CDS', 'intron', 'three_prime_UTR', 'five_prime_UTR', 
                    'nc_RNA', 'miRNA', 'other', 'exon', 'pseudogene', 'lnc_RNA',
                    'pseudogenic_transcript', 'gene', 'mRNA', 'biological_region',
                    'ncRNA_gene', 'scRNA', 'transcript', 'ncRNA_gene-biological_region',
                    'gene_wo_miRNA', 'non_coding_gene', 'protein_coding_gene',
                    'suspected_CDS', 'suspected_exon', 'suspected_five_prime_UTR',
                    'suspected_three_prime_UTR', 'snoRNA']
features_deleted = ['{}/{}.T.{}.bed'.format(bam_folder, research_id, f) for f in features_deleted]
features_deleted = ' '.join(features_deleted)
os.system('rm {}'.format(features_deleted))

os.system('rm {}/{}.bam'.format(bam_folder, research_id))
os.system('rm {}/{}.sorted.bam'.format(bam_folder, research_id))
os.system('rm {}/{}.sorted.bam.bai'.format(bam_folder, research_id))
os.system('rm {}/{}.sorted.bed'.format(bam_folder, research_id))
os.system('rm {}/{}.T.bed'.format(bam_folder, research_id))
os.system('rm {}/{}.T.filtered.bam'.format(bam_folder, research_id))
print('Deleting intermediate files complete!')

# # Peak calling
# chroms = "chr1;chr2"
# peaks_file_name = '{}/{}.T.pureclip.bed'.format(bam_folder, research_id)
# os.system('pureclip --nt 50 --in {} --bai {}.bai --genome {} --iv "{}" --out {}'.format(
#     bam_sorted, bam_sorted, genome_fasta, chroms, peaks_file_name))

# Peak calling
peaks_file_base = '{}/{}.T.peaks'.format(bam_folder, research_id)
peaks_file_name = peaks_file_base + '.bed'
os.system('/bin/bash call_peaks.sh -e {} -b {} -o {} -s {}'.format(
    clipper_env, bam_sorted, peaks_file_base + '.bed', clipper_species))

# Adding peaks to dataframe
print('Adding peaks to dataframe...')

file_base = '{}/{}.T'.format(bam_folder, research_id)
target_sites_file = '{}.filtered.bed'.format(file_base)

df = pd.read_csv(csv_file_base + '.2.csv', sep='\t')
df.insert(len(df.columns), 'peak', None)
df.insert(len(df.columns), 'peak_score', None)

# Add peaks ID to peaks file
peaks_file = open(peaks_file_name, 'r')
peaks_lines = peaks_file.readlines()
peaks_file.close()
for i in range(len(peaks_lines)):
    peaks_lines[i] = peaks_lines[i].split('\t')
    peaks_lines[i][3] = 'peak' + str(i)
    peaks_lines[i] = '\t'.join(peaks_lines[i])

peaks_file = open(peaks_file_name, 'w')
for line in peaks_lines:
    peaks_file.write(line)
peaks_file.close()

# Save genes GTF annotation for coordinates
gene_ids = set(df['gene_id'])
gene_grep = '|'.join(['ID=' + c for c in gene_ids])
ann_format = genome_ann.split('.')[-1]
genes_ann_file = '{}/{}.genes.{}'.format(bam_folder, research_id, ann_format)
os.system('grep -E "{}" {} > {}'.format(
    gene_grep, genome_ann, genes_ann_file))

# Expand peak to minimum 50 nucleotides
peaks = BedTool(peaks_file_name)
genes = BedTool(genes_ann_file)

peaks_scope = 50
peaks = expand_intervals(peaks, peaks_scope, True)
peaks_file_base = '{}/{}.T.peaks'.format(bam_folder, research_id)
peaks_expanded_bed_file = '{}.min{}.bed'.format(peaks_file_base, peaks_scope)
peaks = peaks.saveas(peaks_expanded_bed_file)

peaks_seq = peaks.sequence(fi=genome_fasta, nameOnly=True, s=True)
peaks_expanded_fasta = '{}/{}.T.peaks.min{}.fasta'.format(fastq_folder, research_id, peaks_scope)
peaks_seq.save_seqs(peaks_expanded_fasta)
os.system("sed -i 's/(.*)//' {}".format(peaks_expanded_fasta))

# Add peaks to dataframe
for index, row in df.iterrows():
    target_range = set(range(row['chrom_start'], row['chrom_end']))
    for peak in peaks:
        peak_range = set(range(peak.start, peak.end))
        if row['chrom_strand'] == peak.strand and \
        row['chromosome'] == peak.chrom and \
        bool(peak_range & target_range):
            df.loc[df.index == index, 'peak'] = peak.fields[3]
            df.loc[df.index == index, 'peak_score'] = peak.fields[4]

df = df[df['peak'].notna()]

# Write to log
now = datetime.now().strftime('%Y-%m-%d %H:%M')
os.system('echo "(((((Chimeras after peak calling|||||{})))))\n{}" >> {}'.format(
    now, len(df), log_file))

# Count and drop duplicates
df.insert(len(df.columns), 'count', 1)
df.drop(['mirna_start', 'mirna_end'], axis='columns', inplace=True)
df = df.groupby(['mirna_id', 'is_mirna_first', 'is_revcomp',
                 'feature', 'target_id', 'gene_id', 'gene_type', 'chromosome',
                 'chrom_start', 'chrom_end', 'chrom_strand', 'peak',
                 'peak_score'], dropna=False, as_index=False).size()
df.to_csv(csv_file_base + '.3.csv', index=False, sep='\t')
print('Adding peaks to dataframe complete!')

# Determining miRNA-target interactions
print('Determining miRNA-target interactions...')
mirnas = SeqIO.to_dict(SeqIO.parse(mirna_file, 'fasta'))

df = pd.read_csv(csv_file_base + '.3.csv', sep='\t')
peaks_seq = SeqIO.to_dict(SeqIO.parse(peaks_expanded_fasta, 'fasta'))
peaks_bed = BedTool(peaks_expanded_bed_file)

df.insert(len(df.columns) - 1, 'mfe_rnahybrid', None)
df.insert(len(df.columns) - 1, 'interaction_pattern_rnahybrid', None)

df.insert(len(df.columns) - 1, 'mfe_intarna', None)
df.insert(len(df.columns) - 1, 'interaction_pattern_intarna', None)

df.insert(len(df.columns) - 1, 'mfe_rnafold', None)
df.insert(len(df.columns) - 1, 'interaction_pattern_rnafold', None)

df.insert(len(df.columns) - 1, 'mfe_unafold', None)
df.insert(len(df.columns) - 1, 'interaction_pattern_unafold', None)

for index, row in df.iterrows():
    # RNAhybrid interaction
    rnahybrid_output = subprocess.check_output(['RNAhybrid',
                                                '-s',
                                                '3utr_human',
                                                str(peaks_seq[row['peak']].seq),
                                                str(mirnas[row['mirna_id']].seq)],
                                                encoding='UTF-8').strip().split('\n')
    
    mfe = float(rnahybrid_output[5].replace('mfe:', '').replace('kcal/mol', '').strip())
    df.loc[df.index == index, 'mfe_rnahybrid'] = mfe
    
    target_match = rnahybrid_output[10][10:-3]
    mirna_match = rnahybrid_output[11][10:-3]
    mirna_mismatch = rnahybrid_output[12][10:-3]
    mirna_three_end = mirna_mismatch.find('-')
    mirna_five_end = mirna_mismatch.rfind('-')
    interaction_pattern = ''
    
    for i in range(len(mirna_match)):
        if mirna_match[i] != ' ':
            interaction_pattern += '('
        elif mirna_match[i] == ' ' and mirna_mismatch[i] != ' ':
            interaction_pattern += '.'
    
    interaction_pattern = interaction_pattern[::-1]
    
    df.loc[df.index == index, 'interaction_pattern_rnahybrid'] = interaction_pattern
    
    # Determine interaction position
    position = int(rnahybrid_output[8].replace('position', '').strip()) - 1 # to 0-based
    mirna_len = len(mirnas[row['mirna_id']])
    
    df.loc[df.index == index, 'chrom_start'] = row['chrom_start'] + position
    df.loc[df.index == index, 'chrom_end'] = row['chrom_start'] + position + mirna_len
    
    # IntaRNA interaction
    intarna_output = subprocess.check_output(['IntaRNA',
                                              str(peaks_seq[row['peak']].seq),
                                              str(mirnas[row['mirna_id']].seq)],
                                              encoding='UTF-8').strip()
    if intarna_output != 'no significant hybridization found':
        intarna_output = intarna_output.split('\n')[3:]

        mirna_match = intarna_output[2]
        mirna_mismatch = intarna_output[3]
        mirna_three_end = mirna_mismatch.find('-')
        mirna_five_end = mirna_mismatch.rfind('-')
        interraction_pattern = ''
        interraction_sequence = ''

        for i in range(mirna_three_end+1, mirna_five_end):
            if mirna_mismatch[i] != ' ':
                interraction_pattern += '.'
                interraction_sequence += '.'
            elif len(mirna_match) > i and mirna_match[i] != ' ':
                interraction_pattern += '('
                interraction_sequence += mirna_match[i]
        interraction_pattern = interraction_pattern[::-1]

        df.loc[df.index == index, 'interaction_pattern_intarna'] = interraction_pattern

        mfe = intarna_output[-1].split('energy: ')[1].split(' kcal/mol')[0]
        df.loc[df.index == index, 'mfe_intarna'] = float(mfe)
    
    # RNAfold interaction
    tmp_file_name = 'tmp.fa' 
    tmp_file = open(tmp_file_name, 'w')
    tmp_content = '>tmp\n{}{}'.format(str(mirnas[row['mirna_id']].seq),
                                     peaks_seq[row['peak']].seq)
    tmp_file.write(tmp_content)
    rnafold_output = subprocess.check_output(['RNAfold', tmp_file_name, '--noPS'],
                                             encoding='UTF-8').split('\n')
    if len(rnafold_output) < 2:
        continue
    rnafold_output = rnafold_output[2]
    mfe = float(re.search(r'-?[0-9]+\.[0-9]+', rnafold_output).group(0))
    df.loc[df.index == index, 'mfe_rnafold'] = mfe
    interraction_pattern = rnafold_output.split(' ')[0]
    interraction_pattern = interraction_pattern[:len(mirnas[row['mirna_id']])]
    df.loc[df.index == index, 'interaction_pattern_rnafold'] = interraction_pattern
    
    # UNAfold interaction
    print('mfold SEQ="{}" 1>/dev/null 2>/dev/null'.format(tmp_file_name))
    os.system('mfold SEQ="{}" 1>/dev/null 2>/dev/null'.format(tmp_file_name))
    ct_file_name = tmp_file_name + '.ct'
    if os.path.isfile(ct_file_name):
        unafold_output = subprocess.check_output(['ct2db', ct_file_name],
                                            encoding='UTF-8').split('\n')
        if len(unafold_output) < 3:
            continue
        unafold_output = unafold_output[0:3:2]
        mfe = float(re.search(r'-?[0-9]+\.[0-9]+', unafold_output[0]).group(0))
        df.loc[df.index == index, 'mfe_unafold'] = mfe
        interraction_pattern = unafold_output[1][:len(mirnas[row['mirna_id']])]
        df.loc[df.index == index, 'interaction_pattern_unafold'] = interraction_pattern
    
tmp_file.close()
os.system('rm {}*'.format(tmp_file_name))
df.to_csv(csv_file_base + '.4.csv', index=False, sep='\t')
print('Determining miRNA-target interactions complete!')

print('Saving peaks files...')
df = pd.read_csv(csv_file_base + '.4.csv', sep='\t')
peaks_bed_df = df.loc[:,['chromosome', 'chrom_start',
                         'chrom_end', 'peak',
                         'peak_score', 'chrom_strand']].drop_duplicates(subset=['peak'], keep='first')
peaks_mirna_bed_file = '{}.mirna.bed'.format(peaks_file_base)
peaks_bed_df.to_csv(peaks_mirna_bed_file, sep='\t', index=False, header=False)
peaks_bed = BedTool(peaks_mirna_bed_file)

# Save 200 nt peaks for target availability
peaks_bed = expand_intervals(peaks_bed, 200, True)
peaks_bed = peaks_bed.saveas('{}.200.bed'.format(peaks_file_base))
peaks_seq = peaks_bed.sequence(fi=genome_fasta, nameOnly=True, s=True)
peaks_fasta = '{}/{}.T.peaks.200.fasta'.format(fastq_folder, research_id)
peaks_seq.save_seqs(peaks_fasta)
os.system("sed -i 's/(.*)//' {}".format(peaks_fasta))

# Save peaks ±100 flanking nt for conservation
peaks_bed = BedTool(peaks_mirna_bed_file)
peaks_bed = expand_intervals(peaks_bed, 1000, False)
peaks_bed = peaks_bed.saveas('{}.500flanks.bed'.format(peaks_file_base))

peaks_fasta = '{}/{}.T.peaks.500flanks.fasta'.format(fastq_folder, research_id)
peaks_seq = peaks_bed.sequence(fi=genome_fasta, nameOnly=True, s=True)
peaks_seq.save_seqs(peaks_fasta)
os.system("sed -i 's/(.*)//' {}".format(peaks_fasta))
print('Saving peaks files complete!')

# Target availability with RNAfold
print('Assessing target availability with RNAfold...')
df = pd.read_csv(csv_file_base + '.4.csv', sep='\t')

peak_fasta = '{}/{}.T.peaks.200.fasta'.format(fastq_folder, research_id)
os.system('RNAfold --infile={} --outfile={}.peaks.fold --noPS'.format(
    peak_fasta, research_id))
os.system('mv {}.peaks.fold {}'.format(research_id, fastq_folder))

fold_file = open('{}/{}.peaks.fold'.format(fastq_folder, research_id), 'r')
fold = fold_file.read().split('>')[1:]
fold_file.close()

fold_dict = {}
for f in fold:
    fold_arr = f.split('\n')
    dot_bracket = fold_arr[2].split(' (')[0].strip()
    match_score = (len(dot_bracket) - dot_bracket.count('.')) / len(dot_bracket)
    match_score = round(match_score, 2)
    min_free_energy = fold_arr[2].split(' (')[-1].strip()[:-1]
    fold_dict[fold_arr[0]] = [match_score, min_free_energy]

df.insert(len(df.columns) - 1, 'target_matching_score', None)
df.insert(len(df.columns) - 1, 'target_mfe', None)

# Add availability information to dataframe
for index, row in df.iterrows():
    fold_scores = fold_dict[row['peak']]
    df.loc[df.index == index, 'target_matching_score'] = fold_scores[0]
    df.loc[df.index == index, 'target_mfe'] = fold_scores[1]

df.to_csv(csv_file_base + '.5.csv', index=False, sep='\t')
print('Assessing target availability with RNAfold complete!')

# Conservation
for c in ['mirna', '500flanks']:
    input_file_name = '{}.{}.bed'.format(peaks_file_base, c)
    os.system('/bin/bash conservation.sh -i {} -c {}'.format(
        input_file_name, phylop_folder))

# Adding conservation to dataframe
print('Adding conservation to dataframe...')

df = pd.read_csv(csv_file_base + '.5.csv', sep='\t')

peaks_500flanks_file = '{}.500flanks.bed.conservation.bed'.format(peaks_file_base)
peaks_mirna_file = '{}.mirna.bed.conservation.bed'.format(peaks_file_base)

peaks_500flanks_cons = parse_coservation(peaks_500flanks_file)
peaks_mirna_cons = parse_coservation(peaks_mirna_file)

df.insert(len(df.columns), 'mean_500flanks_conservation', None)
df.insert(len(df.columns), 'mean_mirna_conservation', None)
df.insert(len(df.columns), 'per_base_conservation_500flanks', None)

for index, row in df.iterrows():
    if '' not in peaks_500flanks_cons[row['peak']]:
        flanks_cons = peaks_500flanks_cons[row['peak']][:101] + peaks_500flanks_cons[row['peak']][-100:]
        df.loc[df.index == index, 'mean_500flanks_conservation'] = statistics.mean(flanks_cons)
        per_base_conservation = [str(c) for c in peaks_500flanks_cons[row['peak']]]
        df.loc[df.index == index, 'per_base_conservation_500flanks'] = '|'.join(per_base_conservation)
    if '' not in peaks_mirna_cons[row['peak']]:
        df.loc[df.index == index, 'mean_mirna_conservation'] = statistics.mean(peaks_mirna_cons[row['peak']])
df.to_csv(csv_file_base + '.6.csv', index=False, sep='\t')
        
print('Adding conservation to dataframe complete!')
