#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import argparse
import openpyxl

parser = argparse.ArgumentParser(description='Somatic Mutation Calling: Please see below options')
parser.add_argument("--output-dir", "-o", dest="outputdir", required=True, help="Specify output directory", nargs=1)
parser.add_argument("--script-dir", "-s", dest="scriptdir", required=True, help="Specify directory to files for running Somatic_mutation_calling.py file", nargs=1)
parser.add_argument("--aligned-sample-dir", "-g", dest = "alignedsampledir", required=True, help="Specify directory to folder containing sorted by chromosome coordinates BAM files (i.e. containing *.bam files)", nargs=1)
args = parser.parse_args()

## Setting directory to scripts folder on server where WGS sample info is located
os.chdir(args.scriptdir[0]+"/Files_for_scripts/")

## Set path to folder with aligned BAM files
path_to_bam_files=args.alignedsampledir[0]

## Extract BAM files from the folder

## Run loop to generate file with code to call somatic SNVs and short INDELs for each sample (includes tumor and matched normal)
main_file=open(path_to_bam_files+"/schedule_jobs_var_preprocessing.sh","w")
for sample in glob.glob(path_to_bam_files+"/*.bam"):
    sampleid=os.path.basename(sample).split(".bam")[0]
    file = open(path_to_bam_files+"/"+ sampleid + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH --time=168:00:00 --cpus-per-task=16 --mem-per-cpu=24100 --job-name=" + path_to_bam_files+"/"+ sampleid + ".sh" +
    " -o " + path_to_bam_files + "/" + sampleid+".out -e " + path_to_bam_files+"/"+sampleid+".err\n")
    file.write("export script_dir="+args.scriptdir[0]+"\n")
    file.write("export OD=" + args.outputdir[0] + "\n")
    file.write("export NTHREADS=32\n")
    file.write("export ref_gen=$script_dir/Files_for_scripts/GRCh38.p13.genome.fa\n")
    file.write("export gatk=$script_dir/Packages/gatk/gatk\n")
    file.write("export picard=$script_dir/Packages/picard/build/libs/picard-2.23.1-2-g3f22b39-SNAPSHOT-all.jar\n\n")
    file.write("if [ ! -e $OD/SNVFiles/" + sampleid +"_gatk ];then\n")
    file.write("    echo 'Generating new directory for the sample'" + sampleid + "\n")
    file.write("    mkdir -p $OD/SNVFiles/" + sampleid + "_gatk\n")
    file.write("    cd $OD/SNVFiles/" + sampleid + "_gatk\n\n")
    file.write("    echo 'Marking the duplicate reads'\n")
    file.write("    $gatk MarkDuplicates -I $OD/SNVFiles/" + sampleid + ".bam ")
    file.write("-O $OD/SNVFiles/" + sampleid + "_MD.bam ")
    file.write("-M $OD/SNVFiles/" + sampleid + "_MD.txt\n\n")
    file.write("    echo 'Base score report generation'\n")
    file.write("    $gatk BaseRecalibrator -I $OD/SNVFiles/" + sampleid + "_MD.bam ")
    file.write("-R $ref_gen  --known-sites $script_dir/Files_for_scripts/All_20180418.vcf.gz ")
    file.write("--known-sites $script_dir/Files_for_scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz ")
    file.write("-O " + sampleid + "_MD_BQ1.table\n\n")
    file.write("    echo 'Apply report info to adjust base quality scores'\n")
    file.write("    $gatk ApplyBQSR -R $ref_gen -I " + sampleid "_MD.bam ")
    file.write("--bqsr-recal-file " + sampleid + "_MD_BQ1.table ")
    file.write("-O " + sampleid + "_MD_BQ.bam\n\n")
    file.write("    echo 'Generate index for base recalibrated BAM file\n")
    file.write("    java -jar $picard BuildBamIndex -I " + sampleid + "_MD_BQ.bam ")
    file.write("-O " + sampleid + "_MD_BQ.bai\n\n")
    file.write("    echo 'Base score report and PDF file generation to see the quality of adjustment scores of bases'\n")
    file.write("    $gatk BaseRecalibrator -I " + sampleid + "_MD_BQ.bam ")
    file.write("-R $ref_gen  --known-sites $script_dir/Files_for_scripts/All_20180418.vcf.gz ")
    file.write("--known-sites $script_dir/Files_for_scripts/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz ")
    file.write("-O " + sampleid + "_MD_BQ2.table\n\n")
    file.write("    $gatk AnalyzeCovariates -before " + sampleid + "_MD_BQ1.table ")
    file.write("-after " + sampleid + "_MD_BQ2.table ")
    file.write("-plots " + sampleid + ".pdf\n\n")
    file.write("rm " + sampleid + "_MD.bam\n")
    file.write("fi\n\n")
    file.close()
    main_file.write("sbatch " + path_to_bam_files+"/"+ sampleid + ".sh\n")
main_file.write("chmod +x " + path_to_bam_files+"/schedule_jobs_var_preprocessing.sh")



