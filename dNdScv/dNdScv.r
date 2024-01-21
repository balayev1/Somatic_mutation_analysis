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

dir.now=getwd()
## Build reference database for GRCh38.p13 genome

setwd(reference.dir) # Set directory to the annotation and genome files

# Build reference database
if (!file.exists(paste(reference.dir,"/mutation_output_refcds.rda", sep=""))){
    buildref(cdsfile="Biomart_human_GRCh38.p13.txt", genomefile="GRCh38.p13.genome.fa", outfile = paste(dir.now, "/mutation_output_refcds.rda", sep=""), onlychrs=c(paste("chr", c(rep(1:22,1),'X','Y'), sep="")), useids = T)
    system(paste("sudo cp -r ", dir.now, "/mutation_output_refcds.rda ", reference.dir, sep="")) # transfer RefCDS to output.dir
}

## Run gene mutation analysis

# Load mutation list file
mutations = read.delim(file, header=T)

# Split mutation list into EBV-only and EBV+KSHV-infected samples
ebv.samples = c("E4-3-R1", "E4-3-R2", "E4-3-R3", "E4-9-R1", "E4-9-R2", "E4-9-R3")
ebv.kshv.samples = c("EK4-11-R1", "EK4-11-R2", "EK4-11-R3", "EK3-11-R1", "EK3-11-R2", "EK3-11-R3", "EK3-13-R1", "EK3-13-R2")

ebv.mutations = mutations[mutations$sample_id %in% ebv.samples,]
nrow(ebv.mutations)
# [1] 55857
ebv.kshv.mutations = mutations[mutations$sample_id %in% ebv.kshv.samples,]
nrow(ebv.kshv.mutations)
# [1] 38587


# Run the analysis
## constrain_wnon_wspl = 0 to remove assumption that rate of both nonsense and splicesite mutations must be estimated together
## refdb: reference database
## max_muts_per_gene_per_sample: maximum number of mutations per gene unlimited to detect kataegis
## max_coding_muts_per_sample: maximum number of mutations per sample unlimited to detect kataegis

ebv.dndsout = dndscv(ebv.mutations, refdb = "mutation_output_refcds.rda", sm = "192r_3w", constrain_wnon_wspl = 0, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, cv = NULL, use_indel_sites = F)

ebv.kshv.dndsout = dndscv(ebv.kshv.mutations, refdb = "mutation_output_refcds.rda", sm = "192r_3w", constrain_wnon_wspl = 0, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, cv = NULL, use_indel_sites = F)

dir.create(paste(dir.now, "/dNdScv_IO", sep="")) # Make temprorary directory
setwd(paste(dir.now, "/dNdScv_IO", sep="")) # Set directory to temporary folder

# Save all genes to txt files
write.table(ebv.dndsout$sel_cv, "EBV.huNSG.Sig_Mut.Genes@WRC.Result.txt", sep = "\t", row.names = FALSE)
write.table(ebv.kshv.dndsout$sel_cv, "EBV+KSHV.huNSG.Sig_Mut.Genes@WRC.Result.txt", sep="\t", row.names = FALSE)

# Save global rate parameters
write.table(ebv.dndsout$globaldnds, "EBV.huNSG.Glob.Param@WRC.Result.txt", sep = "\t", row.names = FALSE)
write.table(ebv.kshv.dndsout$globaldnds, "EBV+KSHV.huNSG.Glob.Param@WRC.Result.txt", sep="\t", row.names = FALSE)

# Save sample-wise mutation list matched to gene
write.table(ebv.dndsout$annotmuts, "EBV.huNSG.Sample.Mut@WRC.Result.txt", sep = "\t", row.names = FALSE)
write.table(ebv.kshv.dndsout$annotmuts, "EBV+KSHV.huNSG.Sample.Mut@WRC.Result.txt", sep="\t", row.names = FALSE)

# Transfer all files to the output directory
system(paste("sudo cp -r ", dir.now, "/dNdScv_IO/. ", wd, sep=""))

# Subset genes with q < 0.05
ebv.qsig = ebv.dndsout$sel_cv[ebv.dndsout$sel_cv$qallsubs_cv < 0.05,]
ebv.kshv.qsig = ebv.kshv.dndsout$sel_cv[ebv.kshv.dndsout$sel_cv$qallsubs_cv < 0.05,]
# Number of significantly mutated genes at WRC sites
nrow(ebv.qsig) # for EBV infections
# [1] 82
nrow(ebv.kshv.qsig) # for EBV+KSHV infections
# [1] 97




