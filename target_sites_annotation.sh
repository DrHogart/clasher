#!/bin/bash
usage() {
    echo ""
    echo "Usage: $0 OPTIONS -a <accessions>" 1>&2
    echo ""
    echo "-a - comma-separated accession numbers of reads indexed with bowtie"
    echo "-i - folder with FASTQ files to align"
    echo "-b - folder where to save BED files"
    echo "-l - log file"
    echo "-g - genome index"
    echo "-f - genome GFF3 annotation"
    echo "-m - miRNA hairpins index basename"
    echo ""
}

while getopts ":a:i:b:l:g:f:m:" opt;
do
    case "${opt}" in
        a)
            a=${OPTARG}
            IFS=',' read -r -a accs <<< "$a"
            ;;
        i)
            i=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        g)
            g=${OPTARG}
            ;;
        f)
            f=${OPTARG}
            ;;
        m)
            m=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

echo "Aligning target sites (T) against miRNA hairpins..."
now=$(date +"%Y-%m-%d %H:%M")
echo "(((((Aligning target sites (T) against miRNA|||||$now)))))" >> $l
for acc in "${accs[@]}"
do
    echo $acc
    s="$i/$acc.T.fastq"
    un="$i/$acc.T.nomirna.fastq"
    bowtie --threads 30 -n 1 --sam --seedlen 7 --un $un $m $s 1>/dev/null 2>>$l
done
echo "Aligning target sites (T) against miRNA hairpins complete!"

echo "Aligning target sites (T) against genome..."
for acc in "${accs[@]}"
do
    echo $acc
    now=$(date +"%Y-%m-%d %H:%M")
    echo "(((((Aligning target sites (T) against genome|||||$now)))))" >> $l
    s="$i/$acc.T.nomirna.fastq"
    bam="$b/$acc.T.bam"
    bowtie --threads 10 --sam -m 1 -n 1 --seedlen 8 --best --strata $g $s 2>>$l | samtools view -b - > $bam
done
echo "Aligning target sites (T) against genome complete!"

echo "Intersect T alignment with genome annnotation..."
for acc in "${accs[@]}"
do
    echo $acc
    bam_sorted="$b/$acc.T.bam"
    bed="$b/$acc.T.bed"
    bedtools intersect -s -bed -f 0.8 -wao -a $bam_sorted -b $f > $bed
done
echo "Intersect T alignment with genome annnotation complete!"

echo "Dividing BED files by features..."
for acc in "${accs[@]}"
do
    echo $acc
    features=("lnc_RNA" "pseudogene" "pseudogenic_transcript" "exon" "gene" "mRNA" "biological_region" "ncRNA_gene" "scRNA" "CDS" "miRNA" "transcript" "three_prime_UTR" "five_prime_UTR" "snoRNA")
    name_base="$b/$acc.T"
    tabs=".*	.*	.*	.*	.*	.*	.*	.*	.*	.*	.*	.*	.*	.*	"
    now=$(date +"%Y-%m-%d %H:%M")
    echo "(((((Dividing T alignment by features|||||$now)))))" >> $l
    
    for feature in "${features[@]}"
    do
        echo -n "$feature: " >> $l
        grep "$tabs$feature	" $name_base.bed | wc -l >> $l
        grep "$tabs$feature	" $name_base.bed > $name_base.$feature.bed
    done

    for feature in "exon" "three_prime_UTR" "CDS" "five_prime_UTR"
    do
        mv $name_base.$feature.bed $name_base.suspected_$feature.bed
    done
    
    # gene_wo_miRNA
    bedtools subtract -s -bed -A -a "$name_base.gene.bed" -b "$name_base.miRNA.bed" > \
        "$name_base.gene_wo_miRNA.bed"
    
    # protein_coding_gene
    grep "protein_coding" "$name_base.gene_wo_miRNA.bed" > "$name_base.protein_coding_gene.bed"
    
    # non_coding_gene
    bedtools subtract -s -bed -A -a "$name_base.ncRNA_gene.bed" -b "$name_base.protein_coding_gene.bed" | \
        bedtools subtract -s -bed -A -a - -b "$name_base.miRNA.bed" > "$name_base.non_coding_gene.bed"
    
    # exon
    bedtools intersect -s -bed -u -a "$name_base.suspected_exon.bed" -b "$name_base.protein_coding_gene.bed" > \
        "$name_base.exon.bed"
    
    # three_prime_UTR
    bedtools intersect -s -bed -u -a "$name_base.suspected_three_prime_UTR.bed" -b "$name_base.exon.bed" > \
        "$name_base.three_prime_UTR.bed"
    
    # CDS
    bedtools intersect -s -bed -u -a "$name_base.suspected_CDS.bed" -b "$name_base.exon.bed" | \
        bedtools subtract -s -bed -A -a - -b "$name_base.three_prime_UTR.bed" > "$name_base.CDS.bed"
    
    # five_prime_UTR
    bedtools intersect -s -bed -u -a "$name_base.suspected_five_prime_UTR.bed" -b "$name_base.exon.bed" | \
        bedtools subtract -s -bed -A -a - -b "$name_base.CDS.bed" | \
            bedtools subtract -s -bed -A -a - -b "$name_base.three_prime_UTR.bed" > \
                "$name_base.five_prime_UTR.bed"
                
    # intron
    bedtools subtract -s -bed -A -a "$name_base.protein_coding_gene.bed" -b "$name_base.exon.bed" > \
        "$name_base.intron.bed"
    
    # nc_RNA
    bedtools intersect -s -bed -u -a "$name_base.non_coding_gene.bed" -b "$name_base.lnc_RNA.bed" > "$name_base.nc_RNA.bed"
    bedtools intersect -s -bed -u -a "$name_base.non_coding_gene.bed" -b "$name_base.pseudogene.bed" >> "$name_base.nc_RNA.bed"
    
    # other
    cat "$name_base.ncRNA_gene.bed" > "$name_base.ncRNA_gene-biological_region.bed"
    cat "$name_base.biological_region.bed" >> "$name_base.ncRNA_gene-biological_region.bed"
    bedtools subtract -s -bed -a "$name_base.ncRNA_gene-biological_region.bed" -b "$name_base.protein_coding_gene.bed" | \
    bedtools subtract -s -bed -a - -b "$name_base.non_coding_gene.bed" > "$name_base.other.bed"
done
echo "Dividing BED files by features complete!"