#!/bin/bash
usage() {
    echo ""
    echo "Usage: $0 OPTIONS" 1>&2
    echo ""
    echo "-f - 5'-linker sequence"
    echo "-r - 3'-linker sequence"
    echo "-o - folder for FASTQ/FASTA files"
    echo "-i - folder with index files"
    echo "-b - folder where to save BAM files"
    echo "-c - rRNA index prefix"
    echo "-a - comma-separated accession numbers"
    echo "-m - FASTA file with back transcribed miRNAs"
    echo "-w - research ID"
    echo "-l - log file"
    echo ""
}

while getopts ":f:r:o:i:b:c:a:m:w:l:" opt;
do
    case "${opt}" in
        f)
            f=${OPTARG}
            ;;
        r)
            r=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        i)
            i=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        c)
            c=${OPTARG}
            ;;
        a)
            a=${OPTARG}
            IFS=',' read -r -a accs <<< "$a"
            ;;
        m)
            m=${OPTARG}
            ;;
        w)
            w=${OPTARG}
            ;;
        l)
            l=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

echo "Trimming with Cutadapt..."
now=$(date +"%Y-%m-%d %H:%M")
echo "(((((Trimming with Cutadapt|||||$now)))))" >> $l
arr_xargs=()
for acc in "${accs[@]}"
do
    fastq="$o/$acc.fastq"
    output="$o/$acc.trimmed.fastq"
    arr_xargs+=("cutadapt $fastq --report minimal -g $f -a $r -O 6 --trim-n -m 30 -j 6 --discard-untrimmed -o $output 1>>$l 2>>$l")
done
printf '%s\n' "${arr_xargs[@]}" | xargs --max-procs=20 -i -t sh -c "{}"
echo "Trimming with Cutadapt complete!"

echo "Deduplicating with Prinseq..."
now=$(date +"%Y-%m-%d %H:%M")
echo "(((((Deduplicating with Prinseq|||||$now)))))" >> $l
arr_xargs=()
for acc in "${accs[@]}"
do
    echo $acc
    fastq="$o/$acc.trimmed.fastq"
    output="$o/$acc.dedup"
    prinseq-lite.pl -derep 123 -fastq $fastq -rm_header -out_good $output -out_bad null 2>>$l
done
echo "Deduplicating with Prinseq complete!"

echo "Aligning against rRNA..."
now=$(date +"%Y-%m-%d %H:%M")
echo "(((((Aligning against rRNA|||||$now)))))" >> $l
arr_xargs=()
for acc in "${accs[@]}"
do
    echo $acc
    s="$o/$acc.dedup.fastq"
    un="$o/$acc.norrna.fastq"
    bowtie --threads 30 --un $un $c $s 1>/dev/null 2>>$l
done
echo "Aligning against rRNA complete!"

echo "Reporting quality..."
log_dir=$(dirname $l)
arr_xargs=()
for acc in "${accs[@]}"
do
    file=$o/$acc
    arr_xargs+=("fastqc --quiet -o $log_dir/fastqc $file.fastq")
    arr_xargs+=("fastqc --quiet -o $log_dir/fastqc $file.trimmed.fastq")
#     arr_xargs+=("fastqc --quiet -o $log_dir/fastqc $file.dedup.fastq")
    arr_xargs+=("fastqc --quiet -o $log_dir/fastqc $file.norrna.fastq")
done
printf '%s\n' "${arr_xargs[@]}" | xargs --max-procs=20 -i -t sh -c "{}"
echo "Reporting quality complete!"

echo "Converting FASTQ (trimmed, no rRNA) to FASTA..."
arr_xargs=()
for acc in "${accs[@]}"
do
    arr_xargs+=("seqtk seq -A $o/$acc.norrna.fastq > $o/$acc.norrna.fasta")
done
printf '%s\n' "${arr_xargs[@]}" | xargs --max-procs=20 -i -t sh -c "{}"
echo "Converting FASTQ (trimmed, no rRNA) to FASTA complete!"

echo "Building bowtie indexes..."
mkdir $i/$w
for acc in "${accs[@]}"
do
    echo $acc
    mkdir $i/$w/$acc
    reference_in="$o/$acc.norrna.fasta"
    ebwt_outfile_base="$i/$w/$acc/$acc"
    bowtie-build --quiet --threads 20 $reference_in $ebwt_outfile_base
done
echo "Building bowtie index complete!"

echo "Aligning miRNA against reads..."
now=$(date +"%Y-%m-%d %H:%M")
echo "(((((Aligning miRNA against reads|||||$now)))))" >> $l
for acc in "${accs[@]}"
do
    echo $acc
    ebwt_outfile_base="$i/$w/$acc/$acc"
    bam_base="$b/$acc"
    bowtie --threads 50 --all --sam -f -n 1 --seedlen 8 -e 35 $ebwt_outfile_base $m 2>>$l | samtools view -b - > "$bam_base.bam"
done
echo "Aligning miRNA against reads complete!"

echo "Converting BAM to BED..."
arr_xargs=()
for acc in "${accs[@]}"
do
    bam_base="$b/$acc"
    arr_xargs+=("bedtools bamtobed -i $bam_base.bam >> $bam_base.bed")
done
printf '%s\n' "${arr_xargs[@]}" | xargs --max-procs=20 -i -t sh -c "{}"
echo "Converting BAM to BED complete!"
