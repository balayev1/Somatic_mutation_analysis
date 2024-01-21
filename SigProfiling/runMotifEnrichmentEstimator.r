#!/bin/bash

################ This script generates scripts to run SigProfilerMatrixGenerator.py script for making counts matrix of 
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
              help="full path to MotifEnrichmentEstimator.py script", metavar="character"),
    make_option(c("-g", "--genome-dir"), type="character", default=NULL, 
              help="full path to genome directory", metavar="character"),
    make_option(c("-f", "--genome-name"), type="character", default=NULL, 
              help="genome filename", metavar="character"),
    make_option(c("-m", "--motifs"), default="TC_GA", 
              help="comma separated list of motifs"),
    make_option(c("-c", "--chrs"), default="chr1", 
              help="comma separated list of chromosomes"),       
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for generated files", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

## Generate the output directory structure
output.parent.dir <- opt$'output-dir'
if(file.exists(output.parent.dir)==F){
    cat("Generating top level directory for output...\n")
    dir.create(output.parent.dir)
}

output.dir <- file.path(output.parent.dir, "sample_out")
if(file.exists(output.dir)==F){
    cat("Generating parent directory for sample output...\n")
    dir.create(output.dir)
}

script.dir <- file.path(output.parent.dir, "scripts")
if(file.exists(script.dir)==F){
    cat("Generating script directory...\n")
    dir.create(script.dir)
}

### create a job file
jobsub.filepath <- file.path(script.dir,"motifenrichment.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

### define motifs
motifs <- strsplit(opt$'motifs', ",")[[1]]

## make blank record for processed motifs outside of the loop
completed.sample.names <- vector()

### loop over motifs and generate script to run its enrichment estimator
for (i in motifs){
    ## specify script directory
    script.filepath <- file.path(script.dir, paste("MotifEnrichmentEstimatorfor", i, ".sh", sep = ""))

    ### specify job commands
    cmd.out <- NULL
    cmd.out <- paste0(cmd.out, "#!/bin/bash\n",  
        "#SBATCH --time=24:00:00 --cpus-per-task=1 --mem=128000M --job-name=", paste0("MotifEnrichmentEstimatorfor", i),
            " -o ", file.path(script.dir, paste0("MotifEnrichmentEstimatorfor", i,".out")), " -e ", file.path(script.dir, paste0("MotifEnrichmentEstimatorfor", i, ".err")), "\n")
    
    ### Input
    cmd.out <- paste0(cmd.out, "python -u ", opt$'sig-dir', " --sample-dir ", opt$'vcf-dir', " --genome-dir ", opt$'genome-dir',
        " --genome-file ", opt$'genome-name', " --output-dir ", output.dir, " --motifs ", i, "\n")

    ### write script to the file for running MotifEnrichmentEstimator.py on it
    cat(cmd.out,file=script.filepath,append=F)

    #### Make sure motif was processed before appending sbatch to jobsub.bat
    if (!i %in% completed.sample.names){
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }

    ## record what's so we can check for 
    completed.sample.names <- append(completed.sample.names, i)
    
}

### make it executable
system(paste("chmod 700 ",jobsub.filepath,sep=""))

#################
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " motifs.\n", sep="")
    