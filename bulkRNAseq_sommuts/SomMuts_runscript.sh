#!/bin/bash

################ This script calls the somatic variants as the VCF file using packages picard and GATK in RNA-Seq data
################ It consists of 3 steps: 1) Preprocesses aligned BAM file for variant calling 2) Calls somatic variants 3) Annotation
################ of somatic variants
################ Input: BAM files
################ Output: Annotated VCF files

## Activate conda environment
module load anaconda3
source activate agshin_env

#### Directory with BAM files
export dir_bam=/home/abalay/scratch/Cell_Lines
#### Set dir for outputting the files
export output_dir=/home/abalay/scratch/Cell_Lines/SomMuts

### Please do not change these variables unless they do not exist
export script_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts

## Create directory for somatic mutation analysis
if [ ! -d $output_dir ];then
   mkdir -p $output_dir
fi

############### Step 1: Preprocess aligned BAM file for variant calling
############### Input: BAM file
############### Output: Processed BAM file (marking duplicates, base quality score recalibration)
Rscript $script_dir/SomMutscalling_scripts/SomMuts_preprocessing.r --parent-dir $dir_bam/processed/sample_out --output-dir $output_dir/SomMuts_preprocessing --ref-genome $script_dir/GRCh38.p13.genome.fa --snp-db $script_dir/All_20180418.vcf.gz --indel-db $script_dir/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
$output_dir/SomMuts_preprocessing/scripts/presom/presom.jobsub.bat


############### Step 2: Somatic variant calling
############### Input: Processed BAM file
############### Output: Unfiltered VCF file
Rscript $script_dir/SomMutscalling_scripts/SomMuts_caller.r --parent-dir $output_dir/SomMuts_preprocessing/sample_out --output-dir $output_dir/SomMuts_calling --ref-genome $script_dir/GRCh38.p13.genome.fa --germ-file $script_dir/af-only-gnomad.hg38.vcf.gz --pon  $script_dir/1000g_pon.hg38.vcf.gz --intervals $script_dir/GRCh38.genome.interval_list --anno-file $script_dir/sommuts_anno_file.txt 
$output_dir/SomMuts_calling/scripts/som/som.jobsub.bat



############### Step 3: Somatic variant filtering
############### Input: Processed BAM file + Unfiltered VCF file
############### Output: Filtered VCF file
Rscript $script_dir/SomMutscalling_scripts/SomMuts_variantfilter.r --parent-dir $output_dir/SomMuts_preprocessing/sample_out --output-dir $output_dir/SomMuts_filtering --ref-genome $script_dir/GRCh38.p13.genome.fa  --vcf-dir $output_dir/SomMuts_calling/sample_out --biallelic-file $script_dir/common_all_biallelic.vcf.gz --anno-file $script_dir/sommuts_anno_file.txt 
$output_dir/SomMuts_filtering/scripts/filtersom/filtersom.jobsub.bat


