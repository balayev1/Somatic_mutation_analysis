#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import argparse
import openpyxl

parser = argparse.ArgumentParser(description='Somatic Mutation Calling: Please see below options')
parser.add_argument("--work-dir", "-w", dest="workdir", required=True, help="Specify home directory where all your software is saved", nargs=1)
parser.add_argument("--script-dir", "-s", dest="scriptdir", required=True, help="Specify directory to files for running Somatic_mutation_calling.py file", nargs=1)
parser.add_argument("--gatk-sample-dir", "-g", dest = "gatksampledir", required=True, help="Specify directory to folder containing subfolders with GATK porcessed samples (i.e. containing *_MD_BQ.bam files)", nargs=1)
args = parser.parse_args()

## Setting directory to scripts folder on server where WGS sample info is located
os.chdir(args.scriptdir[0]+"/Files_for_scripts/")
## Download WGS sample info 
xls = pd.ExcelFile('WGS_DBGAP_Sample_Info.xlsx', engine='openpyxl')
## Make a dataframe out of it and preserve 3 vital columns
df1 = pd.read_excel(xls)
df1 = df1[['SRA-ID','Subject_ID','Tumor (yes/no)']]
## Set path to folder with output of "variant_preprocessing" function
path_to_gatk_files=args.gatksampledir[0]
#path_to_gatk_files="/Volumes/G_Cancer_immunobiology$/AGS_AB/WGS_RAW_DLBCL/SNVFiles"
## Run loop to generate file with code to call somatic SNVs and short INDELs for each sample (includes tumor and matched normal)
main_file=open(path_to_gatk_files+"/schedule_jobs.sh","w")
for subid in df1['Subject_ID'].unique():
    df = df1[df1['Subject_ID']==subid]
    file = open(path_to_gatk_files+"/"+ subid + ".sh", "w")
    file.write("#!/bin/bash\n")
    file.write("#SBATCH --time=168:00:00 --cpus-per-task=16 --mem-per-cpu=24100 --job-name=" + path_to_gatk_files+"/"+ subid + ".sh" +
    " -o "+path_to_gatk_files+"/"+subid+".out -e " + path_to_gatk_files+"/"+subid+".err\n")
    file.write("export script_dir="+args.scriptdir[0]+"\n")
#    file.write("export WD=" + args.workdir[0] + "\n")
    file.write("export NTHREADS=32\n")
    file.write("export ref_gen=$script_dir/Files_for_scripts/GRCh38.p13.genome.fa\n")
    file.write("export gatk=$script_dir/Packages/gatk/gatk\n")
    file.write("export picard=$script_dir/Packages/picard/build/libs/picard-2.23.1-2-g3f22b39-SNAPSHOT-all.jar\n\n")
    file.write("mkdir -p " + path_to_gatk_files + "/" + subid + "\n")
    file.write("cd " + path_to_gatk_files + "/" + subid + "\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_unfiltered_mutect2.vcf ] || [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_unfiltered.bamout.bam ];then\n")
    file.write("    $gatk Mutect2 --java-options '-Xmx2G' -R $ref_gen ")
    for histology_index in range(len(df['Tumor (yes/no)'])):
        if df['Tumor (yes/no)'].iloc[histology_index]=='yes': 
            path_to_file_folder=path_to_gatk_files+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_gatk"
            file.write("-I "+ path_to_file_folder+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_sorted_MD_BQ.bam ")
    for histology_index in range(len(df['Tumor (yes/no)'])):
        if df['Tumor (yes/no)'].iloc[histology_index] == 'no':
            path_to_file_folder=path_to_gatk_files+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_gatk"
            file.write("-I "+ path_to_file_folder+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_sorted_MD_BQ.bam ")
    for histology_index in range(len(df['Tumor (yes/no)'])):
        if df['Tumor (yes/no)'].iloc[histology_index] == 'no':
            file.write("-normal " + df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434 ")
    file.write("--f1r2-tar-gz " + subid +"_f1r2.tar.gz ")
    file.write("--germline-resource $script_dir/Files_for_scripts/af-only-gnomad.hg38.vcf.gz ")
    file.write("--panel-of-normals $script_dir/Files_for_scripts/1000g_pon.hg38.vcf.gz ")
    file.write("--native-pair-hmm-threads $NTHREADS ")
    file.write("--af-of-alleles-not-in-resource " + str(0.0000023) + " ")
    file.write("--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter ")
    file.write("--genotype-germline-sites true ")
    file.write("--bam-output " + subid + "_unfiltered.bamout.bam ")
    file.write("--intervals $script_dir/Files_for_scripts/GRCh38.genome.interval_list ")
    file.write("--interval-padding 100 ")
    file.write("-O " + subid + "_unfiltered_mutect2.vcf\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_tumor.pileups.table ];then\n")
    file.write("    $gatk GetPileupSummaries ")
    for histology_index in range(len(df['Tumor (yes/no)'])):
        if df['Tumor (yes/no)'].iloc[histology_index] == 'yes':
            path_to_file_folder=path_to_gatk_files+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_gatk"
            file.write("-I "+ path_to_file_folder+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_sorted_MD_BQ.bam ")
    file.write("-V $script_dir/Files_for_scripts/common_all_biallelic.vcf.gz -L $script_dir/Files_for_scripts/common_all_biallelic.vcf.gz ")
    file.write("-O " + subid + "_tumor.pileups.table\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_normal.pileups.table ];then\n")
    file.write("    $gatk GetPileupSummaries ")
    for histology_index in range(len(df['Tumor (yes/no)'])):
        if df['Tumor (yes/no)'].iloc[histology_index] == 'no':
            path_to_file_folder=path_to_gatk_files+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_gatk"
            file.write("-I "+ path_to_file_folder+"/"+df['SRA-ID'].iloc[histology_index]+"_dbGaP-28434_sorted_MD_BQ.bam ")
    file.write("-V $script_dir/Files_for_scripts/common_all_biallelic.vcf.gz -L $script_dir/Files_for_scripts/common_all_biallelic.vcf.gz ")
    file.write("-O " + subid + "_normal.pileups.table\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + ".contamination.table ] || [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + ".segments.table ];then\n")
    file.write("    $gatk CalculateContamination -I " + subid + "_tumor.pileups.table ")
    file.write("--matched-normal " + subid + "_normal.pileups.table ")
    file.write("-O " + subid + ".contamination.table " + "--tumor-segmentation " + subid + ".segments.table\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_tumor_artifact_prior.tar.gz ];then\n")
    file.write("    $gatk LearnReadOrientationModel -I " + subid +"_f1r2.tar.gz ")
    file.write("-O " + subid + "_tumor_artifact_prior.tar.gz\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + "_filtered_mutect2.vcf ];then\n")
    file.write("    $gatk FilterMutectCalls -R $ref_gen -V " + subid + "_unfiltered_mutect2.vcf ")
    file.write("--contamination-table " + subid + ".contamination.table ")
    file.write("--ob-priors " + subid + "_tumor_artifact_prior.tar.gz " )
    file.write("--tumor-segmentation " + subid + ".segments.table " )
    file.write("-O " + subid + "_filtered_mutect2.vcf\n")
    file.write("fi\n\n")
    file.write("if [ ! -e " + path_to_gatk_files + "/" + subid + "/" + subid + ".artifact_metrics.txt ];then\n")
    file.write("    $gatk CollectSequencingArtifactMetrics -I " + subid + "_unfiltered.bamout.bam ")
    file.write("-R $ref_gen -O " + subid + ".artifact_metrics.txt --DB_SNP $script_dir/Files_for_scripts/All_20180418.vcf.gz\n")
    file.write("fi\n\n")
    file.close()
    main_file.write("sbatch " + path_to_gatk_files+"/"+ subid + ".sh\n")
main_file.write("chmod +x " + path_to_gatk_files+"/schedule_jobs.sh")





