#!/bin/bash

################ This script quantifies dN/dS ratio for missense, nonsense, essential splice site, 5'UTR and 3'UTR mutations
################ It consists of 1 optional and 2 main steps: 1) load reference database (optional) 2) removing duplicate
################ mutations from vcf files of each sample 3) run genome analysis using dNdScv package

## request small number of resources
srun --mem=32GB --time=6:00:00 --pty --cpus-per-task=4 bash

## Activate conda environment
module load anaconda3
source activate sigprofiling_env

## set output directory
export output_dir=/home/abalay/scratch/Cell_Lines/SomMuts/dNdScv_IO/WRC

## set script directory (contains all scripts required to run the pipeline)
export script_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/dNdScv_scripts

## set directory to reference genome FASTA file
export ref_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/GRCh38.p13.genome.fa

Rscript $script_dir/dNdScv_scripts/dNdScv_exec.r \
--genome $ref_dir
--anno $script_dir/Biomart_human_GRCh38.p13.txt
--buildref $script_dir/buildref.r
--dndscv $script_dir/dNdScv.func.r
--file $script_dir/vcfanno.txt
--outputdir $output_dir