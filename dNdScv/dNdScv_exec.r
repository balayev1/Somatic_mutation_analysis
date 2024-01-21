############# dNdScv_exec.r
############# This R script uses dNdScv R package to quantify dN/dS ratios for missense and nonsense mutations and reveal 
############# mutations at essential splice sites, 5'UTR and 3'UTR regions at the level of individual genes. It consists of 3 
############# major steps: 1) creates list of mutations excluding duplicate mutations per sample 2) creates the RefCDS object 
############# from hg38 genome to estimate expected number of mutations at each trinucleotide site per gene. 
############# Note: 192 trinucleotide sites with central mutated 
############# base. Impact of mutations are also estimated (i.e. synonymous, missense, nonsense, splicesite, 5'UTR, 3'UTR).
############# 3) Using negative binomial regression, maximum likelihood of rate parameters for missense, nonsense, splicesite 
############# 5'UTR and 3'UTR mutations is estimated using the computed background mutation rate for each gene(i.e. wmis, 
############# wnon and wspl). Poisson regression computes maximum likelhood of rate parameters at global level.
############# Inputs: full path to genome fasta file
############# full path to file with genomic information on regions of interest (see dNdScv package tutorial)
############# vcf annotation file in TXT format & output directory
############# Output: 
############# Example: Rscript dNdScv_exec.r --reference.dir REFDIR --file ANNOFILE --output.dir OUTPUTDIR

## Load required packages
require(optparse)
require(devtools) 
require(seqinr)
require(dndscv)

## Adjust the arguments 
option_list = list(
    make_option(c("-g", "--genome"), type="character", default=NULL, 
              help="full path to genome FASTA file", metavar="character"),
    make_option(c("-a", "--anno"), type="character", default=NULL, 
              help="full path to genome annotation file", metavar="character"),
    make_option(c("-b", "--buildref"), type="character", default=NULL, 
              help="full path to buildref.r file", metavar="character"),
    make_option(c("-d", "--dndscv"), type="character", default=NULL, 
              help="full path to dNdScv.func.r file", metavar="character")
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="full path to vcf annotation file", metavar="character"),
    make_option(c("-o", "--outputdir"), type="character", default=NULL, 
              help="output directory", metavar="character"));
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

### set output directory
output.dir <- opt$'outputdir'
output.dir <- "/home/abalay/scratch/Cell_Lines/SomMuts/dNdScv_IO/WRC"
if (!dir.exists(output.dir)){
    cat("Creating output directory ...")
    dir.create(output.dir)
}

########### Build reference database for specified genome

### define path to genome FASTA file
genome <- opt$'genome'
genome <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/GRCh38.p13.genome.fa"
### define path to genome annotation file
genanno <- opt$'anno'
genanno <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/Biomart_human_GRCh38.p13.txt"

### specify directory for storing reference database
refdb.dir <- dirname(genome)

### change working directory
setwd(refdb.dir) 

### load buildref.r file
source(opt$'buildref')
source("/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/dNdScv_scripts/buildref.r")

# Build reference database
if (!file.exists(file.path(refdb.dir,"mutation_output_refcds.rda"))){
    buildref(cdsfile=basename(genanno), genomefile=basename(genome), outfile = file.path(refdb.dir, "mutation_output_refcds.rda"), onlychrs=c(paste("chr", c(rep(1:22,1),'X','Y'), sep="")), useids = T)
}

########### Filter duplicate mutations from vcf files
### create input and output subdirectories in output directory for dNdScv run
if (!dir.exists(file.path(output.dir, "Input"))){
    dir.create(file.path(output.dir, "Input"))
}

if (!dir.exists(file.path(output.dir, "Output"))){
    dir.create(file.path(output.dir, "Output"))
}

### open vcf annotation TXT file with three columns as sample ID, full path to corresponding vcf file and condition 
vcfanno <- opt$'file'
vcfanno <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/vcfanno.txt"
vcfanno <- data.frame(read.table(vcfanno, sep="\t"))
colnames(vcfanno) <- c("sample", "vcf", "condition")

### get unique sample IDs
sample.id <- unique(vcfanno$sample)

### create list to add unique mutations per sample
mut.list <- list()
for (i in 1:length(sample.id)){

    cat("Processing vcf files of sample ", sample.id[i], "\n")
    ### get full paths of all vcf files for sample
    samplevcf <- vcfanno$vcf[vcfanno$sample == sample.id[i]]

    sample.list <- list()

    for (j in 1:length(samplevcf)){
        sample.list[[j]] <- read.table(samplevcf, sep="\t", header=FALSE)[,c(1:2,4:5)]
    }

    ### add all mutations to single dataframe
    big.df <- do.call(rbind, sample.list)

    ### check if duplicate mutations exist remove duplicate 
    if (any(duplicated(big.df)) == TRUE) {
        big.df <- big.df %>% distinct()
    }

    ### add all sample mutation to list
    big.df$sample <- rep(sample.id[i], nrow(big.df))
    mut.list[[i]] <- big.df
}

### put all sample mutations into single dataframe
mut.df <- do.call(rbind, mut.list)
mut.df <- mut.df[,c(5,1,2,3,4)]
colnames(mut.df) <- c("sample", "chr", "pos", "ref", "alt")

### save the mutation file
setwd(file.path(output.dir, "Input"))
write.table(mut.df, file="dNdScv_input_mutation_list.txt", sep="\t", row.names=FALSE)

cat("Finished processing duplicate mutations in all samples\n")


########### Run gene mutation analysis
## if there is more than 1 condition, split mutation dataframe 
if (length(unique(vcfanno$condition)) > 1){
    conds <- unique(vcfanno$condition)
    conds.list <- list()
    for (i in 1:length(conds)){
        sample.id <- vcfanno$sample[vcfanno$condition == conds[i]]
        conds.list[[i]] <- mut.df[mut.df$sample %in% sample.id, ] 
        cat("Samples with condition ", conds[i], "contains ", nrow(conds.list[[i]]), "mutations\n")
        names(conds.list)[i] <- conds[i]
    }
}

# Run the analysis
### change working directory 
setwd(refdb.dir) 

### load dNdScv.func.r file
source(opt$'dndscv')
source("/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/dNdScv_scripts/dNdScv.func.r")

## constrain_wnon_wspl = 0 to remove assumption that rate of both nonsense and splicesite mutations must be estimated together
## refdb: reference database
## max_muts_per_gene_per_sample: maximum number of mutations per gene unlimited to detect kataegis
## max_coding_muts_per_sample: maximum number of mutations per sample unlimited to detect kataegis
dndscv.list <- list()
if (length(unique(vcfanno$condition)) > 1){
    for (i in 1:length(conds.list)){
        dndscv.list[[i]] <- dndscv(conds.list[[i]], 
            refdb = "mutation_output_refcds.rda", 
            sm = "192r_3w", 
            constrain_wnon_wspl = 0, 
            max_muts_per_gene_per_sample = Inf, 
            max_coding_muts_per_sample = Inf, 
            cv = NULL, 
            use_indel_sites = F)
        names(dndscv.list)[i] <- names(conds.list)[i]
    }
}
if (length(unique(vcfanno$condition)) == 1){
    dndscv.list[[1]] <- dndscv(mut.df, 
            refdb = "mutation_output_refcds.rda", 
            sm = "192r_3w", 
            constrain_wnon_wspl = 0, 
            max_muts_per_gene_per_sample = Inf, 
            max_coding_muts_per_sample = Inf, 
            cv = NULL, 
            use_indel_sites = F)
}

### change directory
setwd(file.path(output.dir, "Output"))

### Save the output files in output directory
for (i in 1:length(dndscv.list)){
    # Save all genes
    write.table(dndscv.list[[i]]$sel_cv, file=paste0("Sig_Mut.Genes.", names(dndscv.list)[i], ".txt"), sep = "\t", row.names = FALSE)

    # Save global rate parameters
    write.table(dndscv.list[[i]]$globaldnds, file=paste0("Glob.Param.", names(dndscv.list)[i], ".txt"), sep = "\t", row.names = FALSE)

    # Save sample-wise mutation dataframe matched to gene
    write.table(dndscv.list[[i]]$annotmuts, file=paste0("Sample.Mut.", names(dndscv.list)[i], ".txt"), sep = "\t", row.names = FALSE)
}


