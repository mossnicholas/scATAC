#!/bin/bash
module load bedtools
module load R
#SBATCH --mem=10gb
while getopts f:c:o: flag
do
    case "${flag}" in
	f) frags_dir=${OPTARG};;
	c) chip_dir=${OPTARG};;
	o) output_dir=${OPTARG};;
    esac
done

# create mononuc bed files
cd $output_dir
mkdir -p mono_files/
sbatch --mem=30 --wrap="Rscript --vanilla ~/single_nuc/script_test/filter_mono.R $frags_dir $output_dir"
ls mono_files/ > cellIDs.txt
number_cells=$(wc -l < 'cellIDs.txt')

# single-cell coverage across tf sites
mkdir -p ${output_dir}/coverage
for chip_file in $chip_dir/*; do
    chip_name=${chip_file##*/}
        sbatch --array ${number_cells}%100 ~/single_nuc/script_test/nuc_site_coverage.sh ${chip_file} ${chip_name} ${output_dir}
done

 
