#!/bin/bash
usage() {
    echo ""
    echo "Usage: $0 -i <input_bed> -c <chroms_folder>" 1>&2
    echo ""
    echo "<input_bed> - input BED file to be analyzed"
    echo "<chroms_folder> - folder with chromosomes in BED format (converted phyloP widFix files with wig2bed)"
    echo ""
    exit
}

while getopts ":i:c:" opt;
do
    case "${opt}" in
        i)
            i=${OPTARG}
            ;;
        c)
            c=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done


echo "Obtaining conservation..."
chroms=$(awk '{print $1}' $i | sort | uniq)
bed_folder=$(dirname $i)
bed_sorted="$i.sorted.bed"
sort-bed $i > $bed_sorted
bed_per_base="$i.per_base.bed"

awk ' \
    { \
         regionChromosome = $1; \
         regionStart = $2; \
         regionStop = $3; \
         regionID = $4; \
         baseIdx = 0; \
         for (baseStart = regionStart; baseStart < regionStop; baseStart++) { \
             baseStop = baseStart + 1; \
             print regionChromosome"\t"baseStart"\t"baseStop"\t"regionID"-"baseIdx; \
             baseIdx++; \
         } \
    }' $bed_sorted > $bed_per_base
    
arr_cons=()
bed_outputs=()
for chrom in ${chroms[@]}
do
    map_fn="$c/$chrom.bed"
    bed_output="$bed_folder.per_base.$chrom.answer.bed"
    arr_cons+=("bedmap --chrom $chrom --echo --echo-map-score $bed_per_base $map_fn > $bed_output")
    bed_outputs+=($bed_output)
done
printf '%s\n' "${arr_cons[@]}" | xargs --max-procs=24 -i -t sh -c "{}"

cat ${bed_outputs[@]} > "$i.conservation.bed"
rm ${bed_outputs[@]}
rm $bed_sorted
rm $bed_per_base

echo "Obtaining conservation complete!"
echo ""
