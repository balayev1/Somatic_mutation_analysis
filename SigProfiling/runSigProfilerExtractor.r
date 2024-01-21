#!/bin/bash

################ This script generates scripts to run SigProfilerExtractorScript.py script for extracting somatic variant signatures 
################ using non-negative matrix factorization and decomposing to known COSMIC v3.3 signatures
################ 
## Number of cores to multithread
num.cores <- 8

### install required libraries
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-t", "--input-type"), type="character", default=NULL, 
              help="Format of input file(s): 'vcf', 'matrix'", metavar="character"),
    make_option(c("-i", "--input-dir"), type="character", default=NULL, 
              help="Path to directory with VCF files or full path to tab separated file", metavar="character"),
    make_option(c("-s", "--sig-dir"), type="character", default=NULL, 
              help="full path to SigProfilerExtractor.py script", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for generated files", metavar="character"),
    make_option(c("-m", "--max-sig"), type="integer", default=NULL, 
              help="Number of maximum signatures to extract", metavar="integer")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### create output directory if it does not exist
if (!file.exists(opt$'output-dir')){
    cat("Creating output directory ...\n")
    dir.create(opt$'output-dir')
}

### create script directory if it does not exist
script.dir <- file.path(opt$'output-dir', "scripts")
if (!file.exists(script.dir)){
    cat("Creating script directory ...\n")
    dir.create(script.dir)
}

### create file for running SigProfilerExtractor.py
sigextr.filepath<-file.path(script.dir, "sigextr.sh")
if(file.exists(sigextr.filepath)){
    cat("Creating new signature extractor script...\n")
    unlink(sigextr.filepath)
}
file.create(sigextr.filepath)

### specify job commands
cmd.out <- NULL
cmd.out <- paste0(cmd.out, "#!/bin/bash\n",  
    "#SBATCH --time=48:00:00 --cpus-per-task=16 --mem=100000M --job-name=", file.path(script.dir, "SigExtractor.sh"),
        " -o ", file.path(script.dir, "SigExtractor.out"), " -e ", file.path(script.dir, "SigExtractor.err"), "\n")

### Input
cmd.out <- paste0(cmd.out, "python ", opt$'sig-dir', " --input-type ", opt$'input-type', " --output-dir ", opt$'output-dir', 
    " --input-dir ", opt$'input-dir', " --max-signatures ", opt$'max-sig')

### write script to the file for running SigProfilerExtractor.py
cat(cmd.out,file=sigextr.filepath,append=F)

### create a job file
jobsub.filepath <- file.path(script.dir,"sigextr.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

### make it executable
system(paste("chmod 700 ",jobsub.filepath,sep=""))

### append sbatch to job file
cat("sbatch ",sigextr.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
