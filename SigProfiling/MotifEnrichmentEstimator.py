#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################# This script estimates the enrichment levels of variants at specific nucleotide motifs characteristic to members
################# of APOBEC3 family of enzymes (i.e. APOBEC3A-H, AID). As a first step, reference genome is indexed. Then, enrichment levels + 
################# p values for each motif per sample VCF are estimated. In addition, mutations at the selected motifs are extracted
################# and written to a separate vcf file for each sample. 
################# Input: VCF files with variant list. Output: Text file with enrichment level estimates + p_values
################# and VCF files with list of variants at chosen motifs per sample. 

### load required packages
import glob
import os
import numpy as np
import pandas as pd
import re
import itertools
import time
import argparse
from scipy import stats

parser = argparse.ArgumentParser(description='Estimating enrichment levels of mutated nucleotide motifs: Please see below options')
parser.add_argument("--genome-dir", "-g", dest="genomedir", required=True, help="Specify reference genome file directory", nargs=1)
parser.add_argument("--sample-dir", "-s", dest="sampledir", required=True, help="Specify directory to VCF files", nargs=1)
parser.add_argument("--genome-file", "-f", dest = "genomefile", required=True, help="Specify reference genome file name)", nargs=1)
parser.add_argument("--output-dir", "-o", dest = "outputdir", required=True, help="Specify output directory for files)", nargs=1)
parser.add_argument("--chrs", "-c", dest = "chrs", help="Specify list of chromosomes to include for motif enrichment analysis", nargs="+", default=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"])
parser.add_argument("--motifs", "-m", dest = "motifs", type=str, default=['TC_GA','TCW_WGA', 'TCA_TGA', 'YTCW_WGAR','YTCA_TGAR', 'RTCW_WGAY', 'RTCA_TGAY', 'CC_GG','WRC_GYW', 'WRCY_RGYW'], help="Specify list of motifs for enrichment estimation)", nargs="+")
args = parser.parse_args()


####################################################

### Set the variables for genome indexing
GD= args.genomedir[0]
filename = args.genomefile[0]

### Index the genome
def genome_index(filename,GD):	
    os.chdir(GD)
    fasta=open(filename, "r").read().split("\n")
    position_list=[]
    for index in range(len(fasta)):
        if ">" in fasta[index] and "chr" in fasta[index]:
            position_list.append(index)
    chr_dict={}
    for index in range(len(position_list)):
        if fasta[position_list[index]].split(" ")[0][1:] not in chr_dict.keys():
            if index != len(position_list)-1:
                curr_position, next_position = position_list[index], position_list[index+1]
                chr_dict[fasta[curr_position].split(" ")[0][1:]] = "".join(fasta[(curr_position+1):next_position])
            if index==len(position_list)-1:
                curr_position = position_list[index]
                chr_dict[fasta[curr_position].split(" ")[0][1:]]="".join(fasta[(curr_position+1):])
    return chr_dict

start=time.time()
chr_dict=genome_index(filename,GD)
end=time.time()

print(f'Finished indexing the genome in {end-start} second(s)')


#########################################################


### Esimate enrichment level of the needed motifs
def enrichment(vcfdir, output, chr_values, motifs):
    """
    Produces 2 files: 
    1) The file with number of variants, enrichment level and P-Value scores from Fischer's Exact Test'
    2) The VCF file(s) with variants specific to the motif including all the annotation
    Parameters
    _________________
    vcfdir: string
        Path to the vcf files
    output: string
        Output directory for output files
    motifs: list
        List of motifs to be analyzed
        Default=['TC_GA','TCW_WGA', 'TCA_TGA', 'YTCW_WGAR','YTCA_TGAR', 'RTCW_WGAY', 'RTCA_TGAY', 'CC_GG','WRC_GYW', 'WRCY_RGYW']
    genome_index: dictionary
        Dictionary: keys - chromosome ID (e.g. chr1, chr2), values: - genomic sequence
    Output
    _________________
    1) File in 'txt' format
    2) File in 'vcf' format
    """
    ### create output directory if it does not exist
    if os.path.exists(output) == False:
            os.mkdir(output)
    ### start estimation of enrichment of variants at motifs
    if 'TC_GA' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, []   
        ### add variables as keys of dictionary: 
        ### Mut@TC_GA: # of variants at TC_GA motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@TC_GA: # of available TC_GA motifs within ±20 bp vicinity to each TC_GA variant
        ### Enrichment@TC_GA = (Mut@TC_GA/Cont@TC_GA)/(Mut@C_G/Cont@C_G)
        ### P-value_TC_GA: Fischer exact test p-value
        variables_keys=['Mut@TC_GA', 'Mut@C_G', 'Cont@C_G','Cont@TC_GA', 'Enrichment@TC_GA', 'P-value_TC_GA']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at TC_GA motifs
        if os.path.exists(output+"/TC_GA") == False:
            os.mkdir(output+"/TC_GA")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at TC_GA motifs
            output_vcf = output+"/TC_GA"+"/"+variables['SRA-ID'][num] +'_variants@TC_GA.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            tc_ga_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    tc_ga_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_tc_ga, mut_c_g, context_tc_ga, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with TC motifs
                positions_TC = [m.end() for m in re.finditer("TC", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TC = set(positions_TC).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @TC
                mut_tc_ga += len(positions_TC)
                ### find number of TC motifs for each position at TC and sum 
                context_TC = [chr_dict[k][m - 21:m + 20].count('TC') for m in positions_TC]
                context_tc_ga += sum(context_TC)
                ### write variants @TC to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_TC))].to_csv(tc_ga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with GA motifs (reverse complement of TC)
                positions_GA = [m.start()+1 for m in re.finditer("GA", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_GA = set(positions_GA).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @GA
                mut_tc_ga += len(positions_GA)
                ### find number of GA motifs for each position at GA and sum 
                context_GA = [chr_dict[k][m - 21:m + 20].count('GA') for m in positions_GA]
                context_tc_ga += sum(context_GA)
                ### write variants @GA to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_GA))].to_csv(tc_ga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            tc_ga_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            TC_GA_enr=(mut_tc_ga/context_tc_ga)/(mut_c_g/context_c_g)
            ratio1=mut_tc_ga/(mut_c_g-mut_tc_ga); ratio2=context_tc_ga/(context_c_g-context_tc_ga)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_tc_ga, mut_c_g-mut_tc_ga],[context_tc_ga, context_c_g-context_tc_ga]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@TC_GA'].append(mut_tc_ga)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@TC_GA'].append(context_tc_ga)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@TC_GA'].append(TC_GA_enr)
            variables['P-value_TC_GA'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/TC_GA" + "/TC_GA_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for TC_GA motifs in {end-start} second(s)')
    if 'TCW_WGA' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, []  
        ### add variables as keys of dictionary: 
        ### Mut@TCW_WGA: # of variants at TCW_WGA motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@TCW_WGA: # of available TCW_WGA motifs within ±20 bp vicinity to each TCW_WGA variant
        ### Enrichment@TCW_WGA = (Mut@TCW_WGA/Cont@TCW_WGA)/(Mut@C_G/Cont@C_G)
        ### P-value_TCW_WGA: Fischer exact test p-value
        variables_keys=['Mut@TCW_WGA', 'Mut@C_G', 'Cont@C_G','Cont@TCW_WGA', 'Enrichment@TCW_WGA','P-value_TCW_WGA']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at TCW_WGA motifs
        if os.path.exists(output+"/TCW_WGA") == False:
            os.mkdir(output+"/TCW_WGA")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at TCW_WGA motifs
            output_vcf = output+"/TCW_WGA"+"/"+variables['SRA-ID'][num] +'_variants@TCW_WGA.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            tcw_wga_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    tcw_wga_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_tcw_wga, mut_c_g, context_tcw_wga, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with TCW motifs
                positions_TCW = [m.end()-1 for m in re.finditer(r'TCT|TCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TCW = set(positions_TCW).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @TCW
                mut_tcw_wga += len(positions_TCW)
                ### find number of TCW motifs for each position at TCW and sum 
                context_TCA = [chr_dict[k][m - 21:m + 20].count('TCA') for m in positions_TCW]
                context_TCT = [chr_dict[k][m - 21:m + 20].count('TCT') for m in positions_TCW]
                context_tcw_wga += sum(np.add(context_TCA, context_TCT))
                ### write variants @TCW to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_TCW))].to_csv(tcw_wga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with WGA motifs (reverse complement of TCW)
                positions_WGA = [m.end()-1 for m in re.finditer(r'TGA|AGA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_WGA = set(positions_WGA).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @WGA
                mut_tcw_wga += len(positions_WGA)
                ### find number of WGA motifs for each position at WGA and sum 
                context_TGA = [chr_dict[k][m - 21:m + 20].count('TGA') for m in positions_WGA]
                context_AGA = [chr_dict[k][m - 21:m + 20].count('AGA') for m in positions_WGA]
                context_tcw_wga += sum(np.add(context_TGA, context_AGA))
                ### write variants @WGA to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_WGA))].to_csv(tcw_wga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            tcw_wga_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            TCW_WGA_enr=(mut_tcw_wga/context_tcw_wga)/(mut_c_g/context_c_g)
            ratio1=mut_tcw_wga/(mut_c_g-mut_tcw_wga); ratio2=context_tcw_wga/(context_c_g-context_tcw_wga)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_tcw_wga, mut_c_g-mut_tcw_wga],[context_tcw_wga, context_c_g-context_tcw_wga]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@TCW_WGA'].append(mut_tcw_wga)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@TCW_WGA'].append(context_tcw_wga)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@TCW_WGA'].append(TCW_WGA_enr)
            variables['P-value_TCW_WGA'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/TCW_WGA" + "/TCW_WGA_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for TCW_WGA motifs in {end-start} second(s)')
    if 'TCA_TGA' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, []  
        ### add variables as keys of dictionary: 
        ### Mut@TCA_TGA: # of variants at TCA_TGA motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@TCA_TGA: # of available TCA_TGA motifs within ±20 bp vicinity to each TCA_TGA variant
        ### Enrichment@TCA_TGA = (Mut@TCA_TGA/Cont@TCA_TGA)/(Mut@C_G/Cont@C_G)
        ### P-value_TCA_TGA: Fischer exact test p-value
        variables_keys=['Mut@TCA_TGA', 'Mut@C_G', 'Cont@C_G','Cont@TCA_TGA', 'Enrichment@TCA_TGA', 'P-value_TCA_TGA']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at TCA_TGA motifs
        if os.path.exists(output+"/TCA_TGA") == False:
            os.mkdir(output+"/TCA_TGA")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at TCA_TGA motifs
            output_vcf = output+"/TCA_TGA"+"/"+variables['SRA-ID'][num] +'_variants@TCA_TGA.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            tca_tga_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    tca_tga_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_tca_tga, mut_c_g, context_tca_tga, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with TCA motifs
                positions_TCA = [m.end()-1 for m in re.finditer(r'TCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TCA = set(positions_TCA).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @TCA
                mut_tca_tga += len(positions_TCA)
                ### find number of TCA motifs for each position at TCA and sum 
                context_TCA = [chr_dict[k][m - 21:m + 20].count('TCA') for m in positions_TCA]
                context_tca_tga += sum(context_TCA)
                ### write variants @TCA to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_TCA))].to_csv(tca_tga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with TGA motifs (reverse complement of TCA)
                positions_TGA = [m.end()-1 for m in re.finditer(r'TGA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TGA = set(positions_TGA).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @TGA
                mut_tca_tga += len(positions_TGA)
                ### find number of TGA motifs for each position at TGA and sum 
                context_TGA = [chr_dict[k][m - 21:m + 20].count('TGA') for m in positions_TGA]
                context_tca_tga += sum(context_TGA)
                ### write variants @TGA to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_TGA))].to_csv(tca_tga_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            tca_tga_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            TCA_TGA_enr=(mut_tca_tga/context_tca_tga)/(mut_c_g/context_c_g)
            ratio1=mut_tca_tga/(mut_c_g-mut_tca_tga); ratio2=context_tca_tga/(context_c_g-context_tca_tga)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_tca_tga, mut_c_g-mut_tca_tga],[context_tca_tga, context_c_g-context_tca_tga]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@TCA_TGA'].append(mut_tca_tga)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@TCA_TGA'].append(context_tca_tga)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@TCA_TGA'].append(TCA_TGA_enr)
            variables['P-value_TCA_TGA'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/TCA_TGA" + "/TCA_TGA_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for TCA_TGA motifs in {end-start} second(s)')
    if 'YTCW_WGAR' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, []  
        ### add variables as keys of dictionary: 
        ### Mut@YTCW_WGAR: # of variants at YTCW_WGAR motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@YTCW_WGAR: # of available YTCW_WGAR motifs within ±20 bp vicinity to each YTCW_WGAR variant
        ### Enrichment@YTCW_WGAR = (Mut@YTCW_WGAR/Cont@YTCW_WGAR)/(Mut@C_G/Cont@C_G)
        ### P-value_YTCW_WGAR: Fischer exact test p-value
        variables_keys=['Mut@YTCW_WGAR', 'Mut@C_G', 'Cont@C_G','Cont@YTCW_WGAR', 'Enrichment@YTCW_WGAR', 'P-value_YTCW_WGAR']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at YTCW_WGAR motifs
        if os.path.exists(output+"/YTCW_WGAR") == False:
            os.mkdir(output+"/YTCW_WGAR")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at YTCW_WGAR motifs
            output_vcf = output+"/YTCW_WGAR"+"/"+variables['SRA-ID'][num] +'_variants@YTCW_WGAR.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            ytcw_wgar_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    ytcw_wgar_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_ytcw_wgar, mut_c_g, context_ytcw_wgar, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with YTCW motifs
                positions_YTCW = [m.end()-1 for m in re.finditer(r'CTCT|TTCT|CTCA|TTCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_YTCW = set(positions_YTCW).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @YTCW
                mut_ytcw_wgar += len(positions_YTCW)
                ### find number of YTCW motifs for each position at YTCW and sum 
                context_CTCA = [chr_dict[k][m - 21:m + 20].count('CTCA') for m in positions_YTCW]
                context_TTCA = [chr_dict[k][m - 21:m + 20].count('TTCA') for m in positions_YTCW]
                context_CTCT = [chr_dict[k][m - 21:m + 20].count('CTCT') for m in positions_YTCW]
                context_TTCT = [chr_dict[k][m - 21:m + 20].count('TTCT') for m in positions_YTCW]
                array = np.array((context_CTCA, context_TTCA, context_CTCT, context_TTCT))
                context_ytcw_wgar += sum(array.sum(axis=0))
                ### write variants @YTCW to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_YTCW))].to_csv(ytcw_wgar_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with WGAR motifs (reverse complement of YTCW)
                positions_WGAR = [m.end()-2 for m in re.finditer(r'TGAA|TGAG|AGAA|AGAG', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_WGAR = set(positions_WGAR).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @WGAR
                mut_ytcw_wgar += len(positions_WGAR)
                ### find number of WGAR motifs for each position at WGAR and sum 
                context_TGAA = [chr_dict[k][m - 21:m + 20].count('TGAA') for m in positions_WGAR]
                context_TGAG = [chr_dict[k][m - 21:m + 20].count('TGAG') for m in positions_WGAR]
                context_AGAA = [chr_dict[k][m - 21:m + 20].count('AGAA') for m in positions_WGAR]
                context_AGAG = [chr_dict[k][m - 21:m + 20].count('AGAG') for m in positions_WGAR]
                array = np.array((context_TGAA, context_TGAG, context_AGAA, context_AGAG))
                context_ytcw_wgar += sum(array.sum(axis=0))
                ### write variants @WGAR to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_WGAR))].to_csv(ytcw_wgar_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            ytcw_wgar_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            YTCW_WGAR_enr=(mut_ytcw_wgar/context_ytcw_wgar)/(mut_c_g/context_c_g)
            ratio1=mut_ytcw_wgar/(mut_c_g-mut_ytcw_wgar); ratio2=context_ytcw_wgar/(context_c_g-context_ytcw_wgar)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_ytcw_wgar, mut_c_g-mut_ytcw_wgar],[context_ytcw_wgar, context_c_g-context_ytcw_wgar]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@YTCW_WGAR'].append(mut_ytcw_wgar)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@YTCW_WGAR'].append(context_ytcw_wgar)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@YTCW_WGAR'].append(YTCW_WGAR_enr)
            variables['P-value_YTCW_WGAR'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/YTCW_WGAR" + "/YTCW_WGAR_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for YTCW_WGAR motifs in {end-start} second(s)')
    if 'YTCA_TGAR' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@YTCA_TGAR: # of variants at YTCA_TGAR motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@YTCA_TGAR: # of available YTCA_TGAR motifs within ±20 bp vicinity to each YTCA_TGAR variant
        ### Enrichment@YTCA_TGAR = (Mut@YTCA_TGAR/Cont@YTCA_TGAR)/(Mut@C_G/Cont@C_G)
        ### P-value_YTCA_TGAR: Fischer exact test p-value
        variables_keys=['Mut@YTCA_TGAR', 'Mut@C_G', 'Cont@C_G','Cont@YTCA_TGAR', 'Enrichment@YTCA_TGAR', 'P-value_YTCA_TGAR']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at YTCA_TGAR motifs
        if os.path.exists(output+"/YTCA_TGAR") == False:
            os.mkdir(output+"/YTCA_TGAR")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at YTCA_TGAR motifs
            output_vcf = output+"/YTCA_TGAR"+"/"+variables['SRA-ID'][num] +'_variants@YTCA_TGAR.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            ytca_tgar_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    ytca_tgar_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_ytca_tgar, mut_c_g, context_ytca_tgar, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with YTCA motifs
                positions_YTCA = [m.end()-1 for m in re.finditer(r'CTCA|TTCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_YTCA = set(positions_YTCA).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @YTCA
                mut_ytca_tgar += len(positions_YTCA)
                ### find number of YTCA motifs for each position at YTCA and sum 
                context_CTCA = [chr_dict[k][m - 21:m + 20].count('CTCA') for m in positions_YTCA]
                context_TTCA = [chr_dict[k][m - 21:m + 20].count('TTCA') for m in positions_YTCA]
                array = np.array((context_CTCA, context_TTCA))
                context_ytca_tgar += sum(array.sum(axis=0))
                ### write variants @YTCA to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_YTCA))].to_csv(ytca_tgar_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with TGAR motifs (reverse complement of YTCA)
                positions_TGAR = [m.end()-2 for m in re.finditer(r'TGAA|TGAG', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TGAR = set(positions_TGAR).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @TGAR
                mut_ytca_tgar += len(positions_TGAR)
                ### find number of TGAR motifs for each position at TGAR and sum 
                context_TGAA = [chr_dict[k][m - 21:m + 20].count('TGAA') for m in positions_TGAR]
                context_TGAG = [chr_dict[k][m - 21:m + 20].count('TGAG') for m in positions_TGAR]
                array = np.array((context_TGAA, context_TGAG))
                context_ytca_tgar += sum(array.sum(axis=0))
                ### write variants @TGAR to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_TGAR))].to_csv(ytca_tgar_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            ytca_tgar_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            YTCA_TGAR_enr=(mut_ytca_tgar/context_ytca_tgar)/(mut_c_g/context_c_g)
            ratio1=mut_ytca_tgar/(mut_c_g-mut_ytca_tgar); ratio2=context_ytca_tgar/(context_c_g-context_ytca_tgar)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_ytca_tgar, mut_c_g-mut_ytca_tgar],[context_ytca_tgar, context_c_g-context_ytca_tgar]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@YTCA_TGAR'].append(mut_ytca_tgar)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@YTCA_TGAR'].append(context_ytca_tgar)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@YTCA_TGAR'].append(YTCA_TGAR_enr)
            variables['P-value_YTCA_TGAR'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/YTCA_TGAR" + "/YTCA_TGAR_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for YTCA_TGAR motifs in {end-start} second(s)')
    if 'RTCW_WGAY' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@RTCW_WGAY: # of variants at RTCW_WGAY motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@RTCW_WGAY: # of available RTCW_WGAY motifs within ±20 bp vicinity to each RTCW_WGAY variant
        ### Enrichment@RTCW_WGAY = (Mut@RTCW_WGAY/Cont@RTCW_WGAY)/(Mut@C_G/Cont@C_G)
        ### P-RTCW_WGAY: Fischer exact test p-value
        variables_keys=['Mut@RTCW_WGAY', 'Mut@C_G', 'Cont@C_G','Cont@RTCW_WGAY', 'Enrichment@RTCW_WGAY','P-RTCW_WGAY']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at RTCW_WGAY motifs
        if os.path.exists(output+"/RTCW_WGAY") == False:
            os.mkdir(output+"/RTCW_WGAY")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at RTCW_WGAY motifs
            output_vcf = output+"/RTCW_WGAY"+"/"+variables['SRA-ID'][num] +'_variants@RTCW_WGAY.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            rtcw_wgay_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    rtcw_wgay_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_rtcw_wgay, mut_c_g, context_rtcw_wgay, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with RTCW motifs
                positions_RTCW = [m.end()-1 for m in re.finditer(r'ATCT|GTCT|ATCA|GTCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_RTCW = set(positions_RTCW).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @RTCW
                mut_rtcw_wgay += len(positions_RTCW)
                ### find number of RTCW motifs for each position at RTCW and sum 
                context_ATCT = [chr_dict[k][m - 21:m + 20].count('ATCT') for m in positions_RTCW]
                context_GTCT = [chr_dict[k][m - 21:m + 20].count('GTCT') for m in positions_RTCW]
                context_ATCA = [chr_dict[k][m - 21:m + 20].count('ATCA') for m in positions_RTCW]
                context_GTCA = [chr_dict[k][m - 21:m + 20].count('GTCA') for m in positions_RTCW]
                array = np.array((context_ATCT, context_GTCT, context_ATCA, context_GTCA))
                context_rtcw_wgay += sum(array.sum(axis=0))
                ### write variants @RTCW to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_RTCW))].to_csv(rtcw_wgay_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with WGAY motifs (reverse complement of RTCW)
                positions_WGAY = [m.end()-2 for m in re.finditer(r'TGAC|TGAT|AGAC|AGAT', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_WGAY = set(positions_WGAY).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @WGAY
                mut_rtcw_wgay += len(positions_WGAY)
                ### find number of WGAY motifs for each position at WGAY and sum 
                context_TGAC = [chr_dict[k][m - 21:m + 20].count('TGAC') for m in positions_WGAY]
                context_TGAT = [chr_dict[k][m - 21:m + 20].count('TGAT') for m in positions_WGAY]
                context_AGAC = [chr_dict[k][m - 21:m + 20].count('AGAC') for m in positions_WGAY]
                context_AGAT = [chr_dict[k][m - 21:m + 20].count('AGAT') for m in positions_WGAY]
                array = np.array((context_TGAC, context_TGAT, context_AGAC, context_AGAT))
                context_rtcw_wgay += sum(array.sum(axis=0))
                ### write variants @WGAY to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_WGAY))].to_csv(rtcw_wgay_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            rtcw_wgay_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            RTCW_WGAY_enr=(mut_rtcw_wgay/context_rtcw_wgay)/(mut_c_g/context_c_g)
            ratio1=mut_rtcw_wgay/(mut_c_g-mut_rtcw_wgay); ratio2=context_rtcw_wgay/(context_c_g-context_rtcw_wgay)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_rtcw_wgay, mut_c_g-mut_rtcw_wgay],[context_rtcw_wgay, context_c_g-context_rtcw_wgay]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@RTCW_WGAY'].append(mut_rtcw_wgay)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@RTCW_WGAY'].append(context_rtcw_wgay)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@RTCW_WGAY'].append(RTCW_WGAY_enr)
            variables['P-RTCW_WGAY'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/RTCW_WGAY" + "/RTCW_WGAY_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for RTCW_WGAY motifs in {end-start} second(s)')
    if 'RTCA_TGAY' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@RTCA_TGAY: # of variants at RTCA_TGAY motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@RTCA_TGAY: # of available RTCA_TGAY motifs within ±20 bp vicinity to each RTCA_TGAY variant
        ### Enrichment@RTCA_TGAY = (Mut@RTCA_TGAY/Cont@RTCA_TGAY)/(Mut@C_G/Cont@C_G)
        ### P-RTCA_TGAY: Fischer exact test p-value
        variables_keys=['Mut@RTCA_TGAY', 'Mut@C_G', 'Cont@C_G','Cont@RTCA_TGAY', 'Enrichment@RTCA_TGAY', 'P-RTCA_TGAY']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at RTCA_TGAY motifs
        if os.path.exists(output+"/RTCA_TGAY") == False:
            os.mkdir(output+"/RTCA_TGAY")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at RTCA_TGAY motifs
            output_vcf = output+"/RTCA_TGAY"+"/"+variables['SRA-ID'][num] +'_variants@RTCA_TGAY.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            rtca_tgay_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    rtca_tgay_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_rtca_tgay, mut_c_g, context_rtca_tgay, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with RTCA motifs
                positions_RTCA = [m.end()-1 for m in re.finditer(r'ATCA|GTCA', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_RTCA = set(positions_RTCA).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @RTCA
                mut_rtca_tgay += len(positions_RTCA)
                ### find number of RTCA motifs for each position at RTCA and sum 
                context_ATCA = [chr_dict[k][m - 21:m + 20].count('ATCA') for m in positions_RTCA]
                context_GTCA = [chr_dict[k][m - 21:m + 20].count('GTCA') for m in positions_RTCA]
                array = np.array((context_ATCA, context_GTCA))
                context_rtca_tgay += sum(array.sum(axis=0))
                ### write variants @RTCA to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_RTCA))].to_csv(rtca_tgay_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with TGAY motifs (reverse complement of RTCA)
                positions_TGAY = [m.end()-2 for m in re.finditer(r'TGAC|TGAT', chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_TGAY = set(positions_TGAY).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @TGAY
                mut_rtca_tgay += len(positions_TGAY)
                ### find number of TGAY motifs for each position at TGAY and sum 
                context_TGAC = [chr_dict[k][m - 21:m + 20].count('TGAC') for m in positions_TGAY]
                context_TGAT = [chr_dict[k][m - 21:m + 20].count('TGAT') for m in positions_TGAY]
                array = np.array((context_TGAC, context_TGAT))
                context_rtca_tgay += sum(array.sum(axis=0))
                ### write variants @TGAY to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_TGAY))].to_csv(rtca_tgay_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            rtca_tgay_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            RTCA_TGAY_enr=(mut_rtca_tgay/context_rtca_tgay)/(mut_c_g/context_c_g)
            ratio1=mut_rtca_tgay/(mut_c_g-mut_rtca_tgay); ratio2=context_rtca_tgay/(context_c_g-context_rtca_tgay)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_rtca_tgay, mut_c_g-mut_rtca_tgay],[context_rtca_tgay, context_c_g-context_rtca_tgay]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@RTCA_TGAY'].append(mut_rtca_tgay)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@RTCA_TGAY'].append(context_rtca_tgay)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@RTCA_TGAY'].append(RTCA_TGAY_enr)
            variables['P-RTCA_TGAY'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/RTCA_TGAY" + "/RTCA_TGAY_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for RTCA_TGAY motifs in {end-start} second(s)')
    if 'CC_GG' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@CC_GG: # of variants at CC_GG motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@CC_GG: # of available CC_GG motifs within ±20 bp vicinity to each CC_GG variant
        ### Enrichment@CC_GG = (Mut@CC_GG/Cont@CC_GG)/(Mut@C_G/Cont@C_G)
        ### P-value_CC_GG: Fischer exact test p-value
        variables_keys=['Mut@CC_GG', 'Mut@C_G', 'Cont@C_G','Cont@CC_GG', 'Enrichment@CC_GG', 'P-value_CC_GG']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at CC_GG motifs
        if os.path.exists(output+"/CC_GG") == False:
            os.mkdir(output+"/CC_GG")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at CC_GG motifs
            output_vcf = output+"/CC_GG"+"/"+variables['SRA-ID'][num] +'_variants@CC_GG.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            cc_gg_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    cc_gg_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_cc_gg, mut_c_g, context_cc_gg, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with CC motifs
                positions_CC = [m.end() for m in re.finditer("CC", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_CC= set(positions_CC).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @CC
                mut_cc_gg += len(positions_CC)
                ### find number of CC motifs for each position at CC and sum 
                context_CC = [chr_dict[k][m - 21:m + 20].count('CC') for m in positions_CC]
                context_cc_gg += sum(context_CC)
                ### write variants @CC to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_CC))].to_csv(cc_gg_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with GG motifs (reverse complement of CC)
                positions_GG = [m.start()+1 for m in re.finditer("GG", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_GG = set(positions_GG).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @GG
                mut_cc_gg += len(positions_GG)
                ### find number of GG motifs for each position at GG and sum 
                context_GG = [chr_dict[k][m - 21:m + 20].count('GG') for m in positions_GG]
                context_cc_gg += sum(context_GG)
                ### write variants @GG to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_GG))].to_csv(cc_gg_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            cc_gg_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            CC_GG_enr=(mut_cc_gg/context_cc_gg)/(mut_c_g/context_c_g)
            ratio1=mut_cc_gg/(mut_c_g-mut_cc_gg); ratio2=context_cc_gg/(context_c_g-context_cc_gg)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_cc_gg, mut_c_g-mut_cc_gg],[context_cc_gg, context_c_g-context_cc_gg]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@CC_GG'].append(mut_cc_gg)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@CC_GG'].append(context_cc_gg)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@CC_GG'].append(CC_GG_enr)
            variables['P-value_CC_GG'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/CC_GG" + "/CC_GG_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for CC_GG motifs in {end-start} second(s)')
    if 'WRC_GYW' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@WRC_GYW: # of variants at WRC_GYW motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@WRC_GYW: # of available WRC_GYW motifs within ±20 bp vicinity to each WRC_GYW variant
        ### Enrichment@WRC_GYW = (Mut@WRC_GYW/Cont@WRC_GYW)/(Mut@C_G/Cont@C_G)
        ### P-value_WRC_GYW: Fischer exact test p-value
        variables_keys=['Mut@WRC_GYW', 'Mut@C_G', 'Cont@C_G','Cont@WRC_GYW', 'Enrichment@WRC_GYW', 'P-value_WRC_GYW']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at WRC_GYW motifs
        if os.path.exists(output+"/WRC_GYW") == False:
            os.mkdir(output+"/WRC_GYW")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at WRC_GYW motifs
            output_vcf = output+"/WRC_GYW"+"/"+variables['SRA-ID'][num] +'_variants@WRC_GYW.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            wrc_gyw_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    wrc_gyw_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_wrc_gyw, mut_c_g, context_wrc_gyw, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with WRC motifs
                positions_WRC = [m.end() for m in re.finditer("TAC|TGC|AAC|AGC", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_WRC= set(positions_WRC).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @WRC
                mut_wrc_gyw += len(positions_WRC)
                ### find number of WRC motifs for each position at WRC and sum 
                context_TAC = [chr_dict[k][m - 21:m + 20].count('TAC') for m in positions_WRC]
                context_TGC = [chr_dict[k][m - 21:m + 20].count('TGC') for m in positions_WRC]
                context_AAC = [chr_dict[k][m - 21:m + 20].count('AAC') for m in positions_WRC]
                context_AGC = [chr_dict[k][m - 21:m + 20].count('AGC') for m in positions_WRC]
                array = np.array((context_TAC, context_TGC, context_AAC, context_AGC))
                context_wrc_gyw += sum(array.sum(axis=0))
                ### write variants @WRC to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_WRC))].to_csv(wrc_gyw_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with GYW motifs (reverse complement of WRC)
                positions_GYW = [m.start()+1 for m in re.finditer("GCT|GCA|GTT|GTA", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_GYW = set(positions_GYW).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @GYW
                mut_wrc_gyw += len(positions_GYW)
                ### find number of GYW motifs for each position at GYW and sum 
                context_GCT = [chr_dict[k][m - 21:m + 20].count('GCT') for m in positions_GYW]
                context_GCA = [chr_dict[k][m - 21:m + 20].count('GCA') for m in positions_GYW]
                context_GTT = [chr_dict[k][m - 21:m + 20].count('GTT') for m in positions_GYW]
                context_GTA = [chr_dict[k][m - 21:m + 20].count('GTA') for m in positions_GYW]
                array = np.array((context_GCT, context_GCA, context_GTT, context_GTA))
                context_wrc_gyw += sum(array.sum(axis=0))
                ### write variants @GYW to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_GYW))].to_csv(wrc_gyw_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            wrc_gyw_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            WRC_GYW_enr=(mut_wrc_gyw/context_wrc_gyw)/(mut_c_g/context_c_g)
            ratio1=mut_wrc_gyw/(mut_c_g-mut_wrc_gyw); ratio2=context_wrc_gyw/(context_c_g-context_wrc_gyw)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_wrc_gyw, mut_c_g-mut_wrc_gyw],[context_wrc_gyw, context_c_g-context_wrc_gyw]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@WRC_GYW'].append(mut_wrc_gyw)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@WRC_GYW'].append(context_wrc_gyw)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@WRC_GYW'].append(WRC_GYW_enr)
            variables['P-value_WRC_GYW'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/WRC_GYW" + "/WRC_GYW_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for WRC_GYW motifs in {end-start} second(s)')
    if 'WRCY_RGYW' in motifs:
        start=time.time()
        ### create dictionary of variables to be saved as enrichment stats
        variables, variables['SRA-ID']={}, [] 
        ### add variables as keys of dictionary: 
        ### Mut@WRCY_RGYW: # of variants at WRCY_RGYW motifs
        ### Mut@C_G: # of variants at C (or G) bases
        ### Cont@C_G: # of available C (or G) bases within ±20 bp vicinity to each C (or G) variant
        ### Cont@WRCY_RGYW: # of available WRCY_RGYW motifs within ±20 bp vicinity to each WRCY_RGYW variant
        ### Enrichment@WRCY_RGYW = (Mut@WRCY_RGYW/Cont@WRCY_RGYW)/(Mut@C_G/Cont@C_G)
        ### P-value_WRCY_RGYW: Fischer exact test p-value
        variables_keys=['Mut@WRCY_RGYW', 'Mut@C_G', 'Cont@C_G','Cont@WRCY_RGYW', 'Enrichment@WRCY_RGYW', 'P-value_WRCY_RGYW']
        for line in variables_keys:
            variables[line]=[]
        ### create directory to write output of enrichment at WRCY_RGYW motifs
        if os.path.exists(output+"/WRCY_RGYW") == False:
            os.mkdir(output+"/WRCY_RGYW")
        ############ Loop over each vcf file and estimate enrichment variables
        for num, file in enumerate(glob.glob(vcfdir + "/*.vcf"), start=0):
            ### get name of the VCF file prior to ".vcf" suffix
            variables['SRA-ID'].append(file.rsplit("/")[-1].split(".vcf")[0])
            ### load each VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    if line.startswith("#CHROM"):
                        vcf_names = [x for x in line.split('\t')]
                        break
                vcf = pd.read_table(vcf_file, delimiter="\t", comment="#", names = vcf_names)
            vcf_file.close()
            ### keep selected chromosomes
            chr_values = chr_values
            vcf = vcf[vcf["#CHROM"].isin(chr_values)]
            ### subset variants with mutation at bases C and G to T/G and A/C respectively
            vcf_c=vcf[(vcf["REF"] == "C")]
            vcf_c=vcf_c[(vcf_c["ALT"] == "T") | (vcf_c["ALT"] == "G")]
            vcf_g=vcf[(vcf["REF"] == "G")]
            vcf_g=vcf_g[(vcf_g["ALT"] == "A") | (vcf_g["ALT"] == "C")]
            ### open new VCF to write variants at WRCY_RGYW motifs
            output_vcf = output+"/WRCY_RGYW"+"/"+variables['SRA-ID'][num] +'_variants@WRCY_RGYW.vcf'
            if os.path.exists(output_vcf):
                os.remove(output_vcf)
            wrcy_rgyw_file = open(output_vcf,'a')
            ### add VCF header to new VCF file
            with open(file, "r") as vcf_file:
                for line in vcf_file:
                    wrcy_rgyw_file.write(line)
                    if line.startswith("#CHROM"):
                        break
            vcf_file.close()
            ### start counting enrichment variables
            mut_wrcy_rgyw, mut_c_g, context_wrcy_rgyw, context_c_g=0, 0, 0, 0
            for k in chr_values:
                ############ Start with WRCY motifs
                positions_WRCY = [m.end()-1 for m in re.finditer("TACC|TACT|TGCC|TGCT|AACC|AACT|AGCC|AGCT", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_WRCY = set(positions_WRCY).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @WRCY
                mut_wrcy_rgyw += len(positions_WRCY)
                ### find number of WRCY motifs for each position at WRCY and sum 
                context_TACC = [chr_dict[k][m - 21:m + 20].count('TACC') for m in positions_WRCY]
                context_TACT = [chr_dict[k][m - 21:m + 20].count('TACT') for m in positions_WRCY]
                context_TGCC = [chr_dict[k][m - 21:m + 20].count('TGCC') for m in positions_WRCY]
                context_TGCT = [chr_dict[k][m - 21:m + 20].count('TGCT') for m in positions_WRCY]
                context_AACC = [chr_dict[k][m - 21:m + 20].count('AACC') for m in positions_WRCY]
                context_AACT = [chr_dict[k][m - 21:m + 20].count('AACT') for m in positions_WRCY]
                context_AGCC = [chr_dict[k][m - 21:m + 20].count('AGCC') for m in positions_WRCY]
                context_AGCT = [chr_dict[k][m - 21:m + 20].count('AGCT') for m in positions_WRCY]
                array = np.array((context_TACC, context_TACT, context_TGCC, context_TGCT, context_AACC, context_AACT,
                    context_AGCC, context_AGCT))
                context_wrcy_rgyw += sum(array.sum(axis=0))
                ### write variants @WRCY to a vcf file
                vcf_c[(vcf_c["#CHROM"] == k) & (vcf_c["POS"].isin(positions_WRCY))].to_csv(wrcy_rgyw_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with C bases
                positions_C = [m.end() for m in re.finditer("C", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_C = set(positions_C).intersection(set(list(vcf_c["POS"][vcf_c["#CHROM"] == k])))
                ### count number of variants @C
                mut_c_g += len(positions_C)
                ### find number of C bases for each position at C and sum 
                context_C = [chr_dict[k][m - 21:m + 20].count('C') for m in positions_C]
                context_c_g += sum(context_C)
                ############ Continue with RGYW motifs (reverse complement of WRCY)
                positions_RGYW = [m.start()+2 for m in re.finditer("AGCT|GGCT|AGCA|GGCA|AGTT|GGTT|AGTA|GGTA", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_RGYW = set(positions_RGYW).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @RGYW
                mut_wrcy_rgyw += len(positions_RGYW)
                ### find number of RGYW motifs for each position at RGYW and sum 
                context_AGCT = [chr_dict[k][m - 21:m + 20].count('AGCT') for m in positions_RGYW]
                context_GGCT = [chr_dict[k][m - 21:m + 20].count('GGCT') for m in positions_RGYW]
                context_AGCA = [chr_dict[k][m - 21:m + 20].count('AGCA') for m in positions_RGYW]
                context_GGCA = [chr_dict[k][m - 21:m + 20].count('GGCA') for m in positions_RGYW]
                context_AGTT = [chr_dict[k][m - 21:m + 20].count('AGTT') for m in positions_RGYW]
                context_GGTT = [chr_dict[k][m - 21:m + 20].count('GGTT') for m in positions_RGYW]
                context_AGTA = [chr_dict[k][m - 21:m + 20].count('AGTA') for m in positions_RGYW]
                context_GGTA = [chr_dict[k][m - 21:m + 20].count('GGTA') for m in positions_RGYW]
                array = np.array((context_AGCT, context_GGCT, context_AGCA, context_GGCA, context_AGTT, context_GGTT,
                    context_AGTA, context_GGTA))
                context_wrcy_rgyw += sum(array.sum(axis=0))
                ### write variants @RGYW to a vcf file
                vcf_g[(vcf_g["#CHROM"] == k) & (vcf_g["POS"].isin(positions_RGYW))].to_csv(wrcy_rgyw_file, sep="\t", mode='a', header=False, index=False)
                ############ Continue with G bases (complement of C)
                positions_G = [m.end() for m in re.finditer("G", chr_dict[k])]
                ### find intersecting positions in vcf vs genome
                positions_G = set(positions_G).intersection(set(list(vcf_g["POS"][vcf_g["#CHROM"] == k])))
                ### count number of variants @G
                mut_c_g += len(positions_G)
                ### find number of G bases for each position at G and sum 
                context_G = [chr_dict[k][m - 21:m + 20].count('G') for m in positions_G]
                context_c_g += sum(context_G)
            wrcy_rgyw_file.close()
            ### estimate enrichment and p-value using Fischer's exact test
            WRCY_RGYW_enr=(mut_wrcy_rgyw/context_wrcy_rgyw)/(mut_c_g/context_c_g)
            ratio1=mut_wrcy_rgyw/(mut_c_g-mut_wrcy_rgyw); ratio2=context_wrcy_rgyw/(context_c_g-context_wrcy_rgyw)
            if ratio1 <= ratio2:
                p_value=1.0
            else:
                table=np.array([[mut_wrcy_rgyw, mut_c_g-mut_wrcy_rgyw],[context_wrcy_rgyw, context_c_g-context_wrcy_rgyw]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
            print("Finished processing ", file, "\n")
            ### add values into variables list
            variables['Mut@WRCY_RGYW'].append(mut_wrcy_rgyw)
            variables['Mut@C_G'].append(mut_c_g)
            variables['Cont@WRCY_RGYW'].append(context_wrcy_rgyw)
            variables['Cont@C_G'].append(context_c_g)
            variables['Enrichment@WRCY_RGYW'].append(WRCY_RGYW_enr)
            variables['P-value_WRCY_RGYW'].append(p_value)
        ### save as text file
        df = pd.DataFrame(variables)
        df.to_csv(output + "/WRCY_RGYW" + "/WRCY_RGYW_stats.txt", sep="\t", mode='a', header=True, index=False)
        end=time.time()
        print(f'Finished computation of enrichment for WRCY_RGYW motifs in {end-start} second(s)')

start= time.perf_counter()
if __name__ == '__main__':
    vcfdir, output, chr_values, motifs = args.sampledir[0], args.outputdir[0], args.chrs, args.motifs
    enrichment(vcfdir, output, chr_values, motifs)
finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} second(s)')