#!/bin/bash

################ This script generates scripts to run SigProfilerMatrixGeneratorScript.py script for making counts matrix of 
################ different types of somatic variants
################ 
## Number of cores to multithread
num.cores <- 8

### install required libraries
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-v", "--vcf-dir"), type="character", default=NULL, 
              help="full path to directory with vcf files", metavar="character"),
    make_option(c("-s", "--sig-dir"), type="character", default=NULL, 
              help="full path to SigProfilerMatrixGenerator.py script", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for generated files", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)


### create output directory if it does not exist
if (!file.exists(opt$'output-dir')){
    cat("Creating output directory ...\n")
    dir.create(opt$'output-dir')
}

### create file for running SigProfilerMatrixGenerator.py
matgen.filepath<-file.path(opt$'output-dir', "matgen.sh")
if(file.exists(matgen.filepath)){
    cat("Creating new matrix generator script...\n")
    unlink(matgen.filepath)
}
file.create(matgen.filepath)

### specify job commands
cmd.out <- NULL
cmd.out <- paste0(cmd.out, "#!/bin/bash\n",  
    "#SBATCH --time=48:00:00 --cpus-per-task=16 --mem=100000M --job-name=", file.path(opt$'output-dir', "SigProfiling.sh"),
        " -o ", file.path(opt$'output-dir', "SigProfiling.out"), " -e ", file.path(opt$'output-dir', "SigProfiling.err"), "\n")
     
### Input
cmd.out <- paste0(cmd.out, "python ", opt$'sig-dir', " --vcf-dir ", opt$'vcf-dir')

### write script to the file for running SigProfilerMatrixGenerator.py
cat(cmd.out,file=matgen.filepath,append=F)


### create a job file
jobsub.filepath <- file.path(opt$'output-dir',"matgen.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

### make it executable
system(paste("chmod 700 ",jobsub.filepath,sep=""))

### append sbatch to job file
cat("sbatch ",matgen.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)

