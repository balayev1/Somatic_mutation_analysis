###############################################
########################### This script generates matrix of counts for different types of somatic variants
########################### Input: Full Path to directory containing VCF files with somatic variants
########################### Output: 

## load required python libraries
import argparse
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

parser = argparse.ArgumentParser(description='Signature Profiling for Somatic Variants: Please see below options')
parser.add_argument("--vcf-dir", "-s", dest="vcfdir", required=True, help="Specify directory to vcf files with somatic variants", nargs=1)
args = parser.parse_args()

## full path to vcf files
vcfdir = args.vcfdir[0]
vcfdir = "/home/abalay/scratch/Cell_Lines/SomMuts/SigProfiling/sample_out"
## install reference genome GRCh38
genInstall.install('GRCh38', rsync=False, bash=True)

## generate matrix of counts for different types of somatic variants
matrices  = matGen.SigProfilerMatrixGeneratorFunc(vcfdir.rsplit("/")[-1], 
    "GRCh38", 
    vcfdir, 
    exome = False,
    chrom_based = False,
    plot = True,
    seqInfo = True)

