#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################# This script estimates the enrichment levels of mutations at specific nucleotide motifs characteristic to members
################# of APOBEC and ADAR family of enzymes. As a first step, reference genome is indexed. Then, enrichment levels + 
################# p values for each motif per sample are estimated. In addition, mutations at the selected motifs are extracted
################# and written to a separate vcf file for each sample. 
################# Input: VCF files with mutation list. Output: Text file with enrichment level estimates + p_values + q_values
################# and VCF files with list of mutations at chosen motifs per sample. 

import glob
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import urllib.request
import gzip
import time
from multiprocessing import Process, Pool
from itertools import product
import argparse
import sys
import subprocess

# If not root user, elevate to root user
if os.geteuid() == 0:
    print("We're root!")
else:
    print("We're not root.")
    subprocess.call(['sudo', 'python3', *sys.argv])
    sys.exit()

parser = argparse.ArgumentParser(description='Estimating enrichment levels of mutated nucleotide motifs: Please see below options')
parser.add_argument("--genome-dir", "-g", dest="genomedir", required=True, help="Specify reference genome file directory", nargs=1)
parser.add_argument("--sample-dir", "-s", dest="sampledir", required=True, help="Specify directory to VCF files", nargs=1)
parser.add_argument("--genome-file", "-f", dest = "genomefile", required=True, help="Specify reference genome file name)", nargs=1)
parser.add_argument("--output-dir", "-o", dest = "output", required=True, help="Specify output directory for files)", nargs=1)
parser.add_argument("--motifs", "-m", dest = "motifs", required=True, type=list, default=['TCW_WGA','CC_GG','WRC_GYW','WA_TW','AG/CT-GG/CC'], help="Specify list of motifs for enrichment estimation)", nargs=1)
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
def enrichment(WD, output, motifs=['TC_WG','TCW_WGA','YTCW_WGAR','RTCW_WGAY','CC_GG','WRC_GYW','WA_TW','AG/CT-GG/CC']):
    """
    Produces 2 files: 
    1) The file with number of mutations, enrichment level and P-Value scores from Fischer's Exact Test'
    2) The file with mutations specific to the motif including all the annotation
    Parameters
    _________________
    WD: string
        Path to the vcf files
    output: string
        Output directory for output files
    motifs: list
        List of motifs to be analyzed
        Default=['TC_WG','TCW_WGA','YTCW_WGAR','RTCW_WGAY','CC_GG','WRC_GYW','WA_TW','AG/CT-GG/CC']
    genome_index: dictionary
        Dictionary: keys - chromosome ID (e.g. chr1, chr2), values: - genomic sequence
    Output
    _________________
    1) File in 'txt' format
    2) File in 'vcf' format
    """
    os.chdir(WD)
    W,R,Y=['A','T'],['A','G'],['C','T']
    wd=glob.glob("*.vcf")
    variables, variables['SRA-ID']={}, []
    for file in wd:
        variables['SRA-ID'].append(file.split("_")[0])
        if "_" not in file:
            variables['SRA-ID'].append(file.split(".")[0])
    if 'TC_WG' in motifs:
        start=time.time()
        if os.path.exists(output+"/TC_WG") == False:
            os.mkdir(output+"/TC_WG")
        variables_keys=['Mut@TC_WG', 'Mut@C_G', 'Cont@C_G','Cont@TC_WG', 'Enrichment@TC_WG',
                            'Min_estimate of TC_WG', 'P-value_TC_WG']
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,apobec3a_b=open(file, "r"), open(output+"/TC_WG"+"/"+variables['SRA-ID'][num] +'_TC_WG_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=apobec3a_b.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=apobec3a_b.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_tc_wg, mut_c_g, context_tc_wg, context_c_g=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if chr_dict[key][int(new[index][1])-2]=='T':
                                _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_tc_wg+=1
                                bases_tc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(22,41))
                                for value in prox_range:
                                    if bases_tc[value-2:value]=='TC':
                                        context_tc_wg+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if chr_dict[key][int(new[index][1])-2] in W and chr_dict[key][int(new[index][1])]=='A':
                                _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_tc_wg+=1
                                bases_tc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(22,41))
                                for value in prox_range:
                                    if bases_tc[value-1]=='G' and bases_tc[value-2] in W:
                                        context_tc_wg+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentTC_WG=(mut_tc_wg*context_c_g)/(context_tc_wg*mut_c_g)
            ratio1=mut_tc_wg/(mut_c_g-mut_tc_wg); ratio2=context_tc_wg/(context_c_g-context_tc_wg)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_TC_WG=0
            else:
                table=np.array([[mut_tc_wg,mut_c_g-mut_tc_wg],[context_tc_wg,context_c_g-context_tc_wg]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_TC_WG= mut_tc_wg*((enrichmentTC_WG-1)/enrichmentTC_WG)
                else:
                    min_estimate_TC_WG=0
            additional=[]
            additional.extend([mut_tc_wg,mut_c_g,context_c_g,
            context_tc_wg,enrichmentTC_WG,min_estimate_TC_WG,p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_TC_WG']), method = 'BH'); variables['Q-value_TC_WG']=list(p_adjust)
        for index in range(len(variables['Q-value_TC_WG'])):
            if variables['Q-value_TC_WG'][index]>0.05:
                variables['Min_estimate of TC_WG'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for TC_WG motifs in {end-start} second(s)')
    if 'TCW_WGA' in motifs:
        start=time.time()
        if os.path.exists(output+"/TCW_WGA") == False:
            os.mkdir(output+"/TCW_WGA")
        variables_keys=['Mut@TCW_WGA', 'Mut@C_G', 'Cont@C_G','Cont@TCW_WGA', 'Enrichment@TCW_WGA',
                            'Min_estimate of TCW_WGA', 'P-value_TCW_WGA']
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,apobec3a_b=open(file, "r"), open(output+"/TCW_WGA"+"/"+variables['SRA-ID'][num] +'_TCW_WGA_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=apobec3a_b.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=apobec3a_b.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_tcw_wga, mut_c_g, context_tcw_wga, context_c_g=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])-2]=='T' and chr_dict[key][int(new[index][1])] in W:
                                    _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_tcw_wga+=1
                                    bases_tcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(23,41))
                                    for value in prox_range:
                                        if bases_tcw[value-2:value]=='TC' and bases_tcw[value] in W:
                                            context_tcw_wga+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])-2] in W and chr_dict[key][int(new[index][1])]=='A':
                                    _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_tcw_wga+=1
                                    bases_tcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(23,41))
                                    for value in prox_range:
                                        if bases_tcw[value-1]+bases_tcw[value]=='GA' and bases_tcw[value-2] in W:
                                            context_tcw_wga+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentTCW_WGA=(mut_tcw_wga*context_c_g)/(context_tcw_wga*mut_c_g)
            ratio1=mut_tcw_wga/(mut_c_g-mut_tcw_wga); ratio2=context_tcw_wga/(context_c_g-context_tcw_wga)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_TCW_WGA=0
            else:
                table=np.array([[mut_tcw_wga,mut_c_g-mut_tcw_wga],[context_tcw_wga,context_c_g-context_tcw_wga]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_TCW_WGA= mut_tcw_wga*((enrichmentTCW_WGA-1)/enrichmentTCW_WGA)
                else:
                    min_estimate_TCW_WGA=0
            additional=[]
            additional.extend([mut_tcw_wga,mut_c_g,context_c_g,
            context_tcw_wga,enrichmentTCW_WGA,min_estimate_TCW_WGA,p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_TCW_WGA']), method = 'BH'); variables['Q-value_TCW_WGA']=list(p_adjust)
        for index in range(len(variables['Q-value_TCW_WGA'])):
            if variables['Q-value_TCW_WGA'][index]>0.05:
                variables['Min_estimate of TCW_WGA'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for TCW_WGA motifs in {end-start} second(s)')
    if 'YTCW_WGAR' in motifs:
        start=time.time()
        if os.path.exists(output+"/YTCW_WGAR") == False:
            os.mkdir(output+"/YTCW_WGAR")
        variables_keys=['Mut@YTCW_WGAR', 'Mut@C_G', 'Cont@C_G','Cont@YTCW_WGAR', 'Enrichment@YTCW_WGAR',
                            'Min_estimate of YTCW_WGAR', 'P-value_YTCW_WGAR']
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,apobec3a_b=open(file, "r"), open(output+"/YTCW_WGAR"+"/"+variables['SRA-ID'][num] +'_YTCW_WGAR_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=apobec3a_b.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=apobec3a_b.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_ytcw_wgar, mut_c_g, context_ytcw_wgar, context_c_g=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                    if chr_dict[key][int(new[index][1])-2]=='T' and chr_dict[key][int(new[index][1])] in W and chr_dict[key][int(new[index][1])-3] in Y:
                                        _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_ytcw_wgar+=1
                                        bases_ytcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,19))+list(range(23,41))
                                        for value in prox_range:
                                            if bases_ytcw[value-2:value]=='TC' and bases_ytcw[value] in W and bases_ytcw[value-3] in Y:
                                                context_ytcw_wgar+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                    if chr_dict[key][int(new[index][1])-2] in W and chr_dict[key][int(new[index][1])]=='A' and chr_dict[key][int(new[index][1])+1] in R:
                                        _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_ytcw_wgar+=1
                                        bases_ytcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,19))+list(range(23,41))
                                        for value in prox_range:
                                            if bases_ytcw[value-2]+bases_ytcw[value-1]=='GA' and bases_ytcw[value-3] in W and bases_ytcw[value] in R:
                                                context_ytcw_wgar+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentYTCW_WGAR=(mut_ytcw_wgar*context_c_g)/(context_ytcw_wgar*mut_c_g)
            ratio1=mut_ytcw_wgar/(mut_c_g-mut_ytcw_wgar); ratio2=context_ytcw_wgar/(context_c_g-context_ytcw_wgar)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_YTCW_WGAR=0
            else:
                table=np.array([[mut_ytcw_wgar,mut_c_g-mut_ytcw_wgar],[context_ytcw_wgar,context_c_g-context_ytcw_wgar]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_YTCW_WGAR= mut_ytcw_wgar*((enrichmentYTCW_WGAR-1)/enrichmentYTCW_WGAR)
                else:
                    min_estimate_YTCW_WGAR=0
            additional=[]
            additional.extend([mut_ytcw_wgar,mut_c_g,context_c_g,
            context_ytcw_wgar,enrichmentYTCW_WGAR,min_estimate_YTCW_WGAR,p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_YTCW_WGAR']), method = 'BH'); variables['Q-value_YTCW_WGAR']=list(p_adjust)
        for index in range(len(variables['Q-value_YTCW_WGAR'])):
            if variables['Q-value_YTCW_WGAR'][index]>0.05:
                variables['Min_estimate of YTCW_WGAR'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for YTCW_WGAR motifs in {end-start} second(s)')
    if 'RTCW_WGAY' in motifs:
        start=time.time()
        if os.path.exists(output+"/RTCW_WGAY") == False:
            os.mkdir(output+"/RTCW_WGAY")
        variables_keys=['Mut@RTCW_WGAY', 'Mut@C_G', 'Cont@C_G','Cont@RTCW_WGAY', 'Enrichment@RTCW_WGAY',
                            'Min_estimate of RTCW_WGAY', 'P-value_RTCW_WGAY']
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,apobec3a_b=open(file, "r"), open(output+"/RTCW_WGAY"+"/"+variables['SRA-ID'][num] +'_RTCW_WGAY_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=apobec3a_b.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=apobec3a_b.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_rtcw_wgay, mut_c_g, context_rtcw_wgay, context_c_g=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                    if chr_dict[key][int(new[index][1])-2]=='T' and chr_dict[key][int(new[index][1])] in W and chr_dict[key][int(new[index][1])-3] in R:
                                        _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_rtcw_wgay+=1
                                        bases_rtcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,19))+list(range(23,41))
                                        for value in prox_range:
                                            if bases_rtcw[value-2:value]=='TC' and bases_rtcw[value] in W and bases_rtcw[value-3] in R:
                                                context_rtcw_wgay+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                    if chr_dict[key][int(new[index][1])-2] in W and chr_dict[key][int(new[index][1])]=='A' and chr_dict[key][int(new[index][1])+1] in Y:
                                        _=apobec3a_b.write('\t'.join(new[index])+'\n'); mut_rtcw_wgay+=1
                                        bases_rtcw=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,19))+list(range(23,41))
                                        for value in prox_range:
                                            if bases_rtcw[value-2]+bases_rtcw[value-1]=='GA' and bases_rtcw[value-3] in W and bases_rtcw[value] in Y:
                                                context_rtcw_wgay+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentRTCW_WGAY=(mut_rtcw_wgay*context_c_g)/(context_rtcw_wgay*mut_c_g)
            ratio1=mut_rtcw_wgay/(mut_c_g-mut_rtcw_wgay); ratio2=context_rtcw_wgay/(context_c_g-context_rtcw_wgay)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_RTCW_WGAY=0
            else:
                table=np.array([[mut_rtcw_wgay,mut_c_g-mut_rtcw_wgay],[context_rtcw_wgay,context_c_g-context_rtcw_wgay]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_YTCW_WGAR= mut_ytcw_wgar*((enrichmentRTCW_WGAY-1)/enrichmentRTCW_WGAY)
                else:
                    min_estimate_YTCW_WGAR=0
            additional=[]
            additional.extend([mut_rtcw_wgay,mut_c_g,context_c_g,
            context_rtcw_wgay,enrichmentRTCW_WGAY,min_estimate_YTCW_WGAR,p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_RTCW_WGAY']), method = 'BH'); variables['Q-value_RTCW_WGAY']=list(p_adjust)
        for index in range(len(variables['Q-value_RTCW_WGAY'])):
            if variables['Q-value_RTCW_WGAY'][index]>0.05:
                variables['Min_estimate of RTCW_WGAY'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for RTCW_WGAY motifs in {end-start} second(s)')
    if 'CC_GG' in motifs:
        start=time.time()
        if os.path.exists(output+"/CC_GG") == False:
            os.mkdir(output+"/CC_GG")
        variables_keys=['Mut@CC_GG', 'Mut@C_G', 'Cont@C_G','Cont@CC_GG', 'Enrichment@CC_GG',
                            'Min_estimate of CC_GG', 'P-value_CC_GG'] 
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,apobec3g=open(file, "r"), open(output+"/CC_GG"+"/"+variables['SRA-ID'][num] +'_CC_GG_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=apobec3g.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=apobec3g.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_cc_gg, mut_c_g, context_cc_gg, context_c_g=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if chr_dict[key][int(new[index][1])-2]=='C':
                                _=apobec3g.write('\t'.join(new[index])+'\n'); mut_cc_gg+=1
                                bases_cc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                for value in prox_range:
                                    if bases_cc[value-1]+bases_cc[value]=='CC':
                                        context_cc_gg+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                            if chr_dict[key][int(new[index][1])]=='G':
                                _=apobec3g.write('\t'.join(new[index])+'\n'); mut_cc_gg+=1
                                bases_cc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                for value in prox_range:
                                    if bases_cc[value-1]+bases_cc[value]=='GG':
                                        context_cc_gg+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentCC_GG=(mut_cc_gg*context_c_g)/(context_cc_gg*mut_c_g)
            ratio1=mut_cc_gg/(mut_c_g-mut_cc_gg); ratio2=context_cc_gg/(context_c_g-context_cc_gg)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_CC_GG=0
            else:
                table=np.array([[mut_cc_gg,mut_c_g-mut_cc_gg],[context_cc_gg,context_c_g-context_cc_gg]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_CC_GG= mut_cc_gg*((enrichmentCC_GG-1)/enrichmentCC_GG)
                else:
                    min_estimate_CC_GG=0
            additional=[]
            additional.extend([mut_cc_gg,mut_c_g,context_c_g,context_cc_gg,enrichmentCC_GG,min_estimate_CC_GG, p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_CC_GG']), method = 'BH'); variables['Q-value_CC_GG']=list(p_adjust)
        for index in range(len(variables['Q-value_CC_GG'])):
            if variables['Q-value_CC_GG'][index]>0.05:
                variables['Min_estimate of CC_GG'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for CC_GG motifs in {end-start} second(s)')
    if 'WRC_GYW' in motifs:
        start=time.time()
        if os.path.exists(output+"/WRC_GYW") == False:
            os.mkdir(output+"/WRC_GYW")
        variables_keys=['Mut@WRC_GYW', 'Mut@WRCY_RGYW', 'Mut@C_G', 'Cont@C_G','Cont@WRC_GYW', 'Cont@WRCY_RGYW', 'Enrichment@WRC_GYW',
                              'Enrichment@WRCY_RGYW', 'Min_estimate of WRC_GYW', 'Min_estimate of WRCY_RGYW', 'P-value_WRC_GYW', 'P-value_WRCY_RGYW']       
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,aid=open(file, "r"), open(output+"/WRC_GYW"+"/"+variables['SRA-ID'][num] +'_WRC_GYW_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "C" or line.split()[3] == "G":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=aid.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=aid.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_wrc_gyw, mut_c_g, context_wrc_gyw, context_c_g, mut_wrcy_rgyw, context_wrcy_rgyw =0, 0, 0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='C' and (new[index][4]=='T' or new[index][4]=='G'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])-3] in W and chr_dict[key][int(new[index][1])-2] in R:
                                    _=aid.write('\t'.join(new[index])+'\n'); mut_wrc_gyw+=1
                                    bases_wrc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(23,41))
                                    for value in prox_range:
                                        if bases_wrc[value-2] in W and bases_wrc[value-1] in R and bases_wrc[value]=='C':
                                            context_wrc_gyw+=1
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions and (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])-3] in W and chr_dict[key][int(new[index][1])-2] in R and chr_dict[key][int(new[index][1])] in Y:
                                    mut_wrcy_rgyw+=1
                                    bases_wrcy=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,20))+list(range(24,41))
                                    for value in prox_range:
                                        if bases_wrcy[value-3] in W and bases_wrcy[value-2] in R and bases_wrcy[value-1]=='C' and bases_wrcy[value] in Y:
                                            context_wrcy_rgyw+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)];bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='C':
                                context_c_g+=1
                    if new[index][3]=='G' and (new[index][4]=='A' or new[index][4]=='C'):
                        mut_c_g+=1
                        if (key + ":" + str(int(new[index][1])-2)) not in chr_positions:
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])+1] in W and chr_dict[key][int(new[index][1])] in Y:
                                    _=aid.write('\t'.join(new[index])+'\n'); mut_wrc_gyw+=1
                                    bases_wrc=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(2,20))+list(range(23,41))
                                    for value in prox_range:
                                        if bases_wrc[value] in W and bases_wrc[value-1] in Y and bases_wrc[value-2]=='G':
                                            context_wrc_gyw+=1
                            if (key + ":" + str(int(new[index][1])-1)) not in chr_positions and (key + ":" + str(int(new[index][1])+1)) not in chr_positions:
                                if chr_dict[key][int(new[index][1])+1] in W and chr_dict[key][int(new[index][1])-2] in R and chr_dict[key][int(new[index][1])] in Y:
                                    mut_wrcy_rgyw+=1
                                    bases_wrcy=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(3,20))+list(range(24,41))
                                    for value in prox_range:
                                        if bases_wrcy[value] in W and bases_wrcy[value-1] in Y and bases_wrcy[value-2]=='G' and bases_wrcy[value-3] in R:
                                            context_wrcy_rgyw+=1
                        bases_c=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_c=bases_c[0:20]+bases_c[21:]
                        for value in range(len(bases_c)):
                            if bases_c[value]=='G':
                                context_c_g+=1
            enrichmentWRC_GYW=(mut_wrc_gyw*context_c_g)/(context_wrc_gyw*mut_c_g)
            enrichmentWRCY_RGYW=(mut_wrcy_rgyw*context_c_g)/(context_wrcy_rgyw*mut_c_g)
            ratio1=mut_wrc_gyw/(mut_c_g-mut_wrc_gyw); ratio2=context_wrc_gyw/(context_c_g-context_wrc_gyw)
            ratio1y=mut_wrcy_rgyw/(mut_c_g-mut_wrcy_rgyw); ratio2y=context_wrcy_rgyw/(context_c_g-context_wrcy_rgyw)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_WRC_GYW=0
            if ratio1 > ratio2:
                table=np.array([[mut_wrc_gyw,mut_c_g-mut_wrc_gyw],[context_wrc_gyw,context_c_g-context_wrc_gyw]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_WRC_GYW= mut_wrc_gyw*((enrichmentWRC_GYW-1)/enrichmentWRC_GYW)
                else:
                    min_estimate_WRC_GYW=0
            if ratio1y <= ratio2y:
                p_value=1.0; min_estimate_WRCY_RGYW=0
            if ratio1y > ratio2y:
                table=np.array([[mut_wrcy_rgyw,mut_c_g-mut_wrcy_rgyw],[context_wrcy_rgyw,context_c_g-context_wrcy_rgyw]])
                odds_ratio_y, p_value_y = stats.fisher_exact(table,alternative="greater")
                if p_value_y <= 0.05:
                    min_estimate_WRCY_RGYW= mut_wrcy_rgyw*((enrichmentWRCY_RGYW-1)/enrichmentWRCY_RGYW)
                else:
                    min_estimate_WRCY_RGYW=0
            additional=[]
            additional.extend([mut_wrc_gyw,mut_wrcy_rgyw, mut_c_g,context_c_g,context_wrc_gyw,context_wrcy_rgyw, enrichmentWRC_GYW,enrichmentWRCY_RGYW, min_estimate_WRC_GYW, min_estimate_WRCY_RGYW, p_value, p_value_y])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_WRC_GYW']), method = 'BH'); variables['Q-value_WRC_GYW']=list(p_adjust)
        for index in range(len(variables['Q-value_WRC_GYW'])):
            if variables['Q-value_WRC_GYW'][index]>0.05:
                variables['Min_estimate of WRC_GYW'][index]=0
        p_adjust_y = rstats.p_adjust(FloatVector(variables['P-value_WRCY_RGYW']), method = 'BH'); variables['Q-value_WRCY_RGYW']=list(p_adjust_y)
        for index in range(len(variables['Q-value_WRCY_RGYW'])):
            if variables['Q-value_WRCY_RGYW'][index]>0.05:
                variables['Min_estimate of WRCY_RGYW'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for WRC(Y)_(R)GYW motifs in {end-start} second(s)')
    if 'AG/CT-GG/CC' in motifs:
        start=time.time()
        if "TABLE1_hg38.txt.gz" not in GD:
            urllib.request.urlretrieve("http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz", filename=GD+"/TABLE1_hg38.txt.gz")
            rediportal=gzip.open(GD+'/TABLE1_hg38.txt.gz','rt')
        if os.path.exists(output+"/AG_CT") == False:
            os.mkdir(output+"/AG_CT")
        variables_keys=['Mut_AG_CT', 'Mut_A_T', 'Cont_A_T', 'Cont_AG_CT', 'Enrichment@AG_CT', 'Min_estimate of AG_CT', 'P-value_AG_CT']
        for line in variables_keys:
            variables[line]=[]
        list_rna_mut, positions=[], {}
        for line in rediportal.readlines():
            list_rna_mut.append(line.split())
        for index in range(1,len(list_rna_mut)):
            if list_rna_mut[index][0] not in positions.keys():
               positions[list_rna_mut[index][0]]=[]
            else:
               a=str(list_rna_mut[index][1]+list_rna_mut[index][2]+list_rna_mut[index][3])
               positions[list_rna_mut[index][0]].append(a)
        for key in positions.keys():
            positions[key] = set(positions[key])
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,adar=open(file, "r"), open(output+"/AG_CT"+"/"+variables['SRA-ID'][num] +'_ADAR_AG_CT_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "A" or line.split()[3] == "T":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=adar.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=adar.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_ag_ct, context_ag_ct, mut_a_t, context_a_t=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='A' and new[index][4]=='G':
                        mut_a_t+=1
                        if key in positions.keys() and str(new[index][1]+new[index][3]+new[index][4]) in positions[key]:
                            if chr_dict[key][int(new[index][1])] not in chr_positions:
                                if chr_dict[key][int(new[index][1])]=='G':
                                    _=adar.write('\t'.join(new[index])+'\n'); mut_ag_ct+=1
                                    bases_ag=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                    for value in prox_range:
                                        if bases_ag[value-1]+bases_ag[value]=='AG':
                                            context_ag_ct+=1
                        bases_a=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_a=bases_a[0:20]+bases_a[21:]
                        for value in range(len(bases_a)):
                            if bases_a[value]=='A':
                                context_a_t+=1
                    if new[index][3]=='T' and new[index][4]=='C':
                        mut_a_t+=1
                        if key in positions.keys() and str(new[index][1]+new[index][3]+new[index][4]) in positions[key]:
                            if chr_dict[key][int(new[index][1])-2] not in chr_positions:
                                if chr_dict[key][int(new[index][1])-2]=='C':
                                    _=adar.write('\t'.join(new[index])+'\n'); mut_ag_ct+=1
                                    bases_ag=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                    for value in prox_range:
                                        if bases_ag[value-1]+bases_ag[value]=='CT':
                                            context_ag_ct+=1
                        bases_a=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_a=bases_a[0:20]+bases_a[21:]
                        for value in range(len(bases_a)):
                            if bases_a[value]=='T':
                                context_a_t+=1
            enrichmentAG_CT=(mut_ag_ct*context_a_t)/(mut_a_t*context_ag_ct)
            ratio1=mut_ag_ct/(mut_a_t-mut_ag_ct); ratio2=context_ag_ct/(context_a_t-context_ag_ct)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_AG_CT=0
            else:
                table=np.array([[mut_ag_ct,mut_a_t-mut_ag_ct],[context_ag_ct,context_a_t-context_ag_ct]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_AG_CT= mut_ag_ct*((enrichmentAG_CT-1)/enrichmentAG_CT)
                else:
                    min_estimate_AG_CT=0
            additional=[]
            additional.extend([mut_ag_ct, mut_a_t, context_a_t, context_ag_ct, enrichmentAG_CT, min_estimate_AG_CT, p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_AG_CT']), method = 'BH'); variables['Q-value_AG_CT']=list(p_adjust)
        for index in range(len(variables['Q-value_AG_CT'])):
            if variables['Q-value_AG_CT'][index]>0.05:
                variables['Min_estimate of AG_CT'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for AG/CT-GG/CC motifs in {end-start} second(s)')
    if 'WA_TW' in motifs:
        start=time.time()
        if os.path.exists(output+"/WA_TW") == False:
            os.mkdir(output+"/WA_TW")
        if 'AG/CT-GG/CC' not in motifs:
            urllib.request.urlretrieve("http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg38.txt.gz", filename=GD+'/TABLE1_hg38.txt.gz')
            rediportal=gzip.open(GD+'/TABLE1_hg38.txt.gz','rt')
            list_rna_mut, positions=[], {}
            for line in rediportal.readlines():
                list_rna_mut.append(line.split())
            for index in range(1,len(list_rna_mut)):
                if list_rna_mut[index][0] not in positions.keys():
                   positions[list_rna_mut[index][0]]=[]
                else:
                   a=str(list_rna_mut[index][1]+list_rna_mut[index][2]+list_rna_mut[index][3])
                   positions[list_rna_mut[index][0]].append(a)
            for key in positions.keys():
                positions[key] = set(positions[key])
        variables_keys=['Mut@WA_TW', 'Mut@A_T', 'Cont@A_T','Cont@WA_TW', 'Enrichment@WA_TW',
                              'Min_estimate of WA_TW', 'P-value_WA_TW']
        for line in variables_keys:
            variables[line]=[]
        for num, file in enumerate(wd, start=0):
            new = []
            fyle,poleta_theta=open(file, "r"), open(output+"/WA_TW"+"/"+variables['SRA-ID'][num] +'_WA_TW_Mut_List.vcf','a')
            for line in fyle.readlines():
                if len(line.split()) > 3:
                    if line.split()[3] == "A" or line.split()[3] == "T":
                        new.append(line.split())
            fyle.close()
            fyle=open(file, "r").read().split("\n")
            for index in range(len(fyle)):
                if "#CHROM" in fyle[index]:
                    _=poleta_theta.write(''.join(fyle[index])+'\n'); fyle=fyle[index+1:]; break
                _=poleta_theta.write(''.join(fyle[index])+('\n'))
            chr_positions = [i[0] + ":" + i[1] for i in new]
            mut_wa_tw, mut_a_t, context_wa_tw, context_a_t=0, 0, 0, 0
            for index in range(len(new)):
                key = new[index][0]
                if key in chr_dict.keys():
                    if new[index][3]=='A' and (new[index][4]=='T' or new[index][4]=='C' or new[index][4]=='G'):
                        mut_a_t+=1
                        if key in positions.keys() and str(new[index][1]+new[index][3]+new[index][4]) not in positions[key]:
                            if chr_dict[key][int(new[index][1])-2] not in chr_positions:
                                if chr_dict[key][int(new[index][1])-2] in W:
                                    _=poleta_theta.write('\t'.join(new[index])+'\n'); mut_wa_tw+=1
                                    bases_wa=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                    for value in prox_range:
                                        if bases_wa[value-1] in W and bases_wa[value]=='A':
                                            context_wa_tw+=1
                        bases_a=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_a=bases_a[0:20]+bases_a[21:]
                        for value in range(len(bases_a)):
                            if bases_a[value]=='A':
                                context_a_t+=1
                    if new[index][3]=='T' and (new[index][4]=='A' or new[index][4]=='G' and new[index][4]=='C'):
                        mut_a_t+=1
                        if key in positions.keys() and str(new[index][1]+new[index][3]+new[index][4]) not in positions[key]:
                            if chr_dict[key][int(new[index][1])] not in chr_positions:
                                if chr_dict[key][int(new[index][1])] in W:
                                    _=poleta_theta.write('\t'.join(new[index])+'\n'); mut_wa_tw+=1
                                    bases_wa=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; prox_range=list(range(1,20))+list(range(22,41))
                                    for value in prox_range:
                                        if bases_wa[value] in W and bases_wa[value-1]=='T':
                                            context_wa_tw+=1
                        bases_a=chr_dict[key][(int(new[index][1])-21):(int(new[index][1])+20)]; bases_a=bases_a[0:20]+bases_a[21:]
                        for value in range(len(bases_a)):
                            if bases_a[value]=='T':
                                context_a_t+=1
            enrichmentWA_TW=(mut_wa_tw*context_a_t)/(context_wa_tw*mut_a_t)
            ratio1=mut_wa_tw/(mut_a_t-mut_wa_tw); ratio2=context_wa_tw/(context_a_t-context_wa_tw)
            if ratio1 <= ratio2:
                p_value=1.0; min_estimate_WA_TW=0
            else:
                table=np.array([[mut_wa_tw,mut_a_t-mut_wa_tw],[context_wa_tw,context_a_t-context_wa_tw]])
                odds_ratio, p_value = stats.fisher_exact(table,alternative="greater")
                if p_value <= 0.05:
                    min_estimate_WA_TW= mut_wa_tw*((enrichmentWA_TW-1)/enrichmentWA_TW)
                else:
                    min_estimate_WA_TW=0
            additional=[]
            additional.extend([mut_wa_tw, mut_a_t, context_a_t, context_wa_tw, enrichmentWA_TW, min_estimate_WA_TW, p_value])
            for index in range(len(variables_keys)):
                variables[variables_keys[index]].append(additional[index])
        rstats = importr('stats')
        p_adjust = rstats.p_adjust(FloatVector(variables['P-value_WA_TW']), method = 'BH'); variables['Q-value_WA_TW']=list(p_adjust)
        for index in range(len(variables['Q-value_WA_TW'])):
            if variables['Q-value_WA_TW'][index]>0.05:
                variables['Min_estimate of WA_TW'][index]=0
        end=time.time()
        print(f'Finished computation of enrichment for WA_TW motifs in {end-start} second(s)')
    df = pd.DataFrame(variables, columns=variables.keys())
    df.to_csv(output+"/"+"Mut.Enrichment.Results.txt", sep="\t", index=False, mode="a")

start= time.perf_counter()

if __name__ == '__main__':
    WD, output, motifs=args.sampledir[0], args.output[0], args.motifs[0]
    enrichment(WD, output)

finish = time.perf_counter()

print(f'Finished in {round(finish-start,2)} second(s)')







