#!/bin/bash
module load bedtools
bed_file="$3/mono_files/$(sed -n ${SLURM_ARRAY_TASK_ID}p '$3/cellIDs.txt')"
bed_name=${bed_file##*/}
bedtools coverage -a $1 -b ${bed_file} -d | awk 'OFS="\t" {print $1":"$2"-"$3, $4, $5}' | awk '{sums[$2] += $3} END { for (i in sums) printf("%s %s\n", i, sums[i])}' | sort -n > $3/coverage/$2_${bed_name}
