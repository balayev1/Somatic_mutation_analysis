###############################################
########################### This script extracts de novo somatic variant signatures from the samples and decomposes to known
########################### COSMIC v3.3 signatures
########################### Input: Full Path to directory containing VCF files or Full Path to tab separated file 
########################### Output: Report on extracted signatures

## load required python libraries
import argparse
import os
from SigProfilerExtractor import sigpro as sig

parser = argparse.ArgumentParser(description='Signature Profiling for Somatic Variants: Please see below options')
parser.add_argument("--input-type", "-t", dest="inputtype", type=str, required=True, help="Specify format of input files (Must be 'vcf' or 'matrix')", nargs=1)
parser.add_argument("--output-dir", "-o", dest="outputdir", type=str, required=True, help="Specify full output directory", nargs=1)
parser.add_argument("--input-dir", "-i", dest="inputdir", type=str, required=True, help="Specify directory if input type is 'vcf' or full path to 'matrix' file", nargs=1)
parser.add_argument("--max-signatures", "-m", dest="maxsig", type=int, required=True, help="Maximum number of signatures to extract", nargs=1)
args = parser.parse_args()

## input type
input_type = args.inputtype[0]

## path to input file(s)
input_path = args.inputdir[0]

## set current working directory
os.chdir(args.outputdir[0])

## set maximum signature as integer
maxsig = int(args.maxsig[0])

## extract signatures from somatic variants
def main_function():
    sig.sigProfilerExtractor(input_type, "sample_out", input_path, minimum_signatures=1, maximum_signatures=maxsig,
        reference_genome = "GRCh38", opportunity_genome="GRCh38", context_type="96", nmf_replicates=500, stability=0.8, 
        cosmic_version= 3.3, make_decomposition_plots=True,
        collapse_to_SBS96=True, exome=False)

if __name__=="__main__":
   main_function()


