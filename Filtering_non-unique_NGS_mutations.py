############ Goal of this script is to create a dataframe with the list of unique mutations from each EBV and EBV+KSHV infected
############ humanized mice (i.e. huNSG mice) sample. In other terms, mutation repeated > 1 in the technical replicates of 
############ each sample is removed. Output: dataframe with list of combined unique mutations from each sample.

## Upload the required packages
import pandas as pd
import os
import glob
import argparse

parser = argparse.ArgumentParser(description='Filtering samplewise non-unique mutations: Please see below options')
parser.add_argument("--annotation-dir", "-a", dest="annotationdir", required=True, help="Specify file directory with sample information Excel file", nargs=1)
parser.add_argument("--sample-dir", "-s", dest="sampledir", required=True, help="Specify directory to VCF files", nargs=1)
parser.add_argument("--annotation-file", "-f", dest = "annfile", required=True, help="Specify annotation file name in Excel format)", nargs=1)
args = parser.parse_args()


## Upload excel file with sample info
os.chdir(args.annotationdir[0]) # Set directory to file with the sample information

print("Uploading sample information file ...")
samples_info = pd.read_excel(args.annfile[0], sheet_name="EBV KSHV in vivo samples") # Upload excel file with the sample information
samples_info = samples_info.iloc[0: , :]  
samples_info.columns = samples_info.iloc[0] 

## Generate dataframe with sample info and respective mutation info
os.chdir(args.sampledir[0]) # Set directory to VCF files

print("Generating dataframe with the mutation list of all samples ...")
list_of_files=glob.glob('*.vcf') # Find any files with ".vcf" filename ending in the VCF directory
vcf_list=[]
for file in list_of_files:
    sra_id=file.split("_")[0] # sample SRA-ID (specific to every replicate)
    index_of_sra_id=samples_info.index[samples_info['SRA-ID']==sra_id].tolist()[0] # sample SRA-ID line index in sample info file
    sample_id=samples_info["Cell line"][index_of_sra_id]+"-"+samples_info["Replicate #"][index_of_sra_id] # Sample-IDs 
    vcf=open(file,'r') # Load sample VCF file
    # Identify the starting line after which mutations are annotated in VCF file
    for num, line in enumerate(vcf,1):
        if "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in line.strip():
            mutation_start_index = num+1 # starting line
    _=vcf.seek(0) # Scroll back to start of VCF file to reread it
    # Extract Sample-ID, chromosome, position, reference and mutant nucleotide(s)
    for line in vcf.readlines()[mutation_start_index:]:
        chrom, pos, ref, alt = line.split("\t")[0], line.split("\t")[1], line.split("\t")[3], line.split("\t")[4]
        line=sample_id+"\t"+chrom+"\t"+pos+"\t"+ref+"\t"+alt
        vcf_list.append(line.split("\t"))
# Dataframe with mutations:
## For each mutation:
### Format: Sample-ID, chromosome #, position on chromosome, reference nucleotide(s) and mutated variant of the reference
vcf_pd=pd.DataFrame(vcf_list, columns=["sample_id", "chrom", "pos", "ref", "mut"])
vcf_pd.shape[0] # Number of mutations
# 112485

## Filter out duplicate mutations in each sample
print("Filtering duplicate mutations in each sample ...")

vcf_list = []
for sample in vcf_pd["sample_id"].unique():
    mini_df = vcf_pd[vcf_pd["sample_id"] == sample] # Subset dataframe by sample ID
    # If duplicate mutations are in the dataframe, keep only one of the duplicates
    if any(mini_df.duplicated()) is True:
        mini_df = mini_df.drop_duplicates(subset=None, keep="first")
    vcf_list.append(mini_df)

# Make dataframe without duplicated mutations in every sample
vcf_pd = pd.concat(vcf_list, axis=0, ignore_index=True)
vcf_pd.shape[0] # Number of mutations (i.e. minus the duplicates)
# 94444



    










    

            




    



