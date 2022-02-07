############# Script 1) creates the RefCDS object from hg38 genome where coding (+splicesite) sequences and expected number of 
############# mutations at each trinucleotide site are estimated per gene. Note: 192 trinucleotide sites with central mutated 
############# base. Impact of mutations are also estimated (i.e. synonymous, missense, nonsense, splicesite).
############# 2) Using negative binomial regression, maximum likelhood of rate parameters for missense, nonsense and splicesite 
############# mutations is estimated using the computed background mutation rate for each gene(i.e. wmis, 
############# wnon and wspl). Poisson regression computes maximum likelhood of rate parameters at global level. 

## Load required packages
library("optparse")
library("dndscv") 
library("seqinr")

## Adjust the arguments 
option_list = list(
  make_option(c("-r", "--reference.dir"), type="character", default=NULL, 
              help="reference genome & annotation file directory", metavar="character"),
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="full path to mutation list file ", metavar="character"),
    make_option(c("-o", "--output.dir"), type="character", default=NULL, 
              help="output directory", metavar="character")
);
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

reference.dir = opt$reference.dir
wd = opt$output.dir
file = opt$file

## Build reference database for GRCh38.p13 genome

setwd(reference.dir) # Set directory to the annotation and genome files

# Build reference database
if (!file.exists(paste(wd,"/mutation_output_refcds.rda", sep=""))){
    buildref(cdsfile="Biomart_human_GRCh38.p13.txt", genomefile="GRCh38.p13.genome.fa", outfile = paste(wd, "/mutation_output_refcds.rda", sep=""), onlychrs=c(paste("chr", c(rep(1:22,1),'X','Y'), sep="")), useids = T)
}

## Run gene mutation analysis

setwd(wd) # Set output directory

# Load mutation list file
mutations = read.delim(file, header=T)

# Run the analysis
## constrain_wnon_wspl = 0 to remove assumption that rate of both nonsense and splicesite mutations must be estimated together
## refdb: reference database
## max_muts_per_gene_per_sample: maximum number of mutations per gene unlimited to detect kataegis
## max_coding_muts_per_sample: maximum number of mutations per sample unlimited to detect kataegis

dndsout = dndscv(mutations, refdb = "mutation_output_refcds.rda", sm = "192r_3w", constrain_wnon_wspl = 0, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, cv = NULL, use_indel_sites = F)


