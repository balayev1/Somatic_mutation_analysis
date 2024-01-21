#!/bin/bash

################ This script runs Signature Profiling for different types of somatic variants
################ As an additional option, enrichment of variants at specific signature motifs can be estimated 
################ (i.e. at TC/WG, TCW/WGA, TCA/TGA, YTCW/WGAR, YTCA/TGAR, RTCW/WGAY, RTCA/TGAY, CC/GG, WRC/GYW, WRCY/RGYW motifs)
################ !!!!!!!!!!!! Note: install SigProfilerMatrixGenerator package (e.g. pip install SigProfilerMatrixGenerator)

## Activate conda environment
module load anaconda3
source activate sigprofiling_env


## set directory to VCF files and output directory
export VCFDIR=/home/abalay/scratch/Cell_Lines/SomMuts/SigProfiling/sample_out
export output_dir=/home/abalay/scratch/Cell_Lines/SomMuts/SigProfiling/scripts

## set script directory (contains all scripts required to run the pipeline)
export script_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/SigProfiling_scripts

############### Step 1: Remove mitochondrial variants from VCF files (Optional)
############### Input: VCF file
############### Output: VCF file without mitochondrial variants
for vcf in $(find $VCFDIR -name "*vcf")
do
    sample=$(basename ${vcf%.*})
    cat $vcf | grep -v '^chrM' > $VCFDIR/"$sample"_nochrM.vcf
    rm $vcf
done

############### Step 2: Generate matrix of counts for different types of somatic variants
############### Input: full path to VCF files + full path to SigProfilerMatrixGeneratorScript.py file
############### Output: Counts Matrix + Plots
Rscript $script_dir/runSigProfilerMatrixGenerator.r --vcf-dir $VCFDIR --sig-dir $script_dir/SigProfilerMatrixGeneratorScript.py --output-dir $output_dir
$output_dir/matgen.jobsub.bat

### if running motif enrichment analysis switch conda environments:
source activate motifenrichment_env

############### Step 3: Run motif enrichment script for each VCF file (Optional)
############### Input: full path to VCF files + full path to reference genome FASTA file + list of motifs + list of chromosomes
############### Output: TXT file per motif + VCF files with variants at motifs
export output_dir_for_motifs=/home/abalay/scratch/Cell_Lines/SomMuts/MotifEnrichment
export genome_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts
Rscript $script_dir/runMotifEnrichmentEstimator.r  --vcf-dir $VCFDIR/input \
--sig-dir $script_dir/MotifEnrichmentEstimator.py \
--output-dir $output_dir_for_motifs \
--genome-dir $genome_dir \
--genome-name GRCh38.p13.genome.fa \
--motifs TC_GA,TCW_WGA,TCA_TGA,YTCW_WGAR,YTCA_TGAR,RTCW_WGAY,RTCA_TGAY,CC_GG,WRC_GYW,WRCY_RGYW
$output_dir_for_motifs/scripts/motifenrichment.jobsub.bat

############### !!! Note: there is also available chromosome specification command which can be envisaged by running command
############### (Rscript $script_dir/runMotifEnrichmentEstimator.sh -h)



############### Step 4: Run de novo signature extraction using somatic variant counts and decompose to COSMIC v3.3 signatures
############### Input: full path to SBS96 matrix from Step 2 output + full path to SigProfilerExtractorScript.py file 
############### Output: Reports on decomposed signatures
############### !!! Note: before running signature extraction, please install SigProfilerExtractor (e.g. pip install SigProfilerExtractor)
export input_file=/home/abalay/scratch/Cell_Lines/SomMuts/SigProfiling/sample_out/output/SBS/sample_out.SBS96.all
export output_dir=/home/abalay/scratch/Cell_Lines/SomMuts/SigProfiling/SigExtraction
Rscript $script_dir/runSigProfilerExtractor.r --input-type matrix \
--sig-dir $script_dir/SigProfilerExtractorScript.py \
--input-dir $input_file \
--output-dir $output_dir \
--max-sig 20
$output_dir/scripts/sigextr.jobsub.bat


############### Step 5: Run enrichment analysis of signature mutations at the level of individual genes using dNdScv package
############### Input: full path to genome FASTA file + full path to genome annotation file + full path to vcf annotation file
############### Output:
############### !!! Note: before running dNdScv, please install
Rscript $script_dir/dNdScv_exec.r
--genome 
--anno
--file

