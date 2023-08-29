#!/bin/bash
usage() {
    echo ""
    echo "Usage: $0 -e <clipper_env> -b <input_bam> -o <output_bed> -s <clipper_species>" 1>&2
    echo ""
    echo "<clipper_env> - conda environment with CLIPper installed"
    echo "<input_bam> - input BAM file"
    echo "<output_bed> - output BED file with peaks"
    echo "<clipper_species> - species; equivalent to `-s` CLIPper option"
    echo ""
    exit
}

while getopts ":e:b:o:s:" opt;
do
    case "${opt}" in
        e)
            e=${OPTARG}
            ;;
        b)
            b=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        s)
            s=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

echo "Peak calling..."
eval "$(conda shell.bash hook)"
conda activate $e
clipper -v -b $b -o $o -s $s
echo "Peak calling complete!"