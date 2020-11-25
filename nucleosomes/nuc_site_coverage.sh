#!/bin/bash
module load bedtools
motif_file=$1
bed_file=$2
outdir=$3
motif_name=${motif_file##*/}
bed_name=${bed_file##*/}

# total coverage across all sites
bedtools coverage -a ${motif_file} -b ${bed_file} -d | awk 'OFS="\t" {print $1":"$2"-"$3, $4, $5}' | awk '{sums[$2] += $3} END { for (i in sums) printf("%s %s\n", i, sums[i])}' | sort -n > ${outdir}/${bed_name}_${motif_name}.txt
