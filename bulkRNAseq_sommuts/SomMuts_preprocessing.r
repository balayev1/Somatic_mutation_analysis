#!/bin/bash

################ This script generates scripts to run preprocessing steps of somatic variant calling

## Number of cores to multithread
num.cores <- 8

### install required libraries
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-p", "--parent-dir"), type="character", default=NULL, 
              help="parent directory with aligned bam files", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for generated output files", metavar="character"),
    make_option(c("-r", "--ref-genome"), type="character", default=NULL, 
              help="path to reference genome", metavar="character"),
    make_option(c("-s", "--snp-db"), type="character", default=NULL, 
              help="path to SNP file ", metavar="character"),
    make_option(c("-d", "--indel-db"), type="character", default=NULL, 
              help="path to germline indel file ", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- opt$p
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}

## Generate the output directory structure
output.parent.dir <- opt$o
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

presom.script.dir <- file.path(script.dir, "presom")
  if(file.exists(presom.script.dir)==F){
      cat("Generating presom script directory...\n")
      dir.create(presom.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(presom.script.dir,"presom.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

### set the path to gatk executable
gatk <- "/home/abalay/scratch/Pipelines_and_Files/Packages/gatk/build/install/gatk/bin/gatk"

### set the path to picard executable
picard <- "/home/abalay/scratch/Pipelines_and_Files/Packages/picard/build/libs/picard-2.23.1-2-g3f22b39-SNAPSHOT-all.jar"

## Extract the paths for the data parent directories (sample folders)
soup.sample.names <- list.files(input.parent.dir)
soup.sample.dir.paths <- file.path(input.parent.dir, soup.sample.names)

## make blank record for processed samples outside of the loop
completed.sample.names <- vector()
## Loop over sample folders
for(i in 1:length(soup.sample.dir.paths)){

    bam.file <- grep(".bam", list.files(soup.sample.dir.paths[i]), value = TRUE)[1]

    sample.output.dir <- file.path(output.dir, soup.sample.names[i])
    if (!dir.exists(sample.output.dir)){
        dir.create(sample.output.dir)
    }

    script.filepath <- file.path(presom.script.dir, paste(soup.sample.names[i], "_presom.sh", sep = ""))

################################### Print the script
    cmd.out <- NULL
    cmd.out <- paste("#!/bin/bash\n\n")
    cmd.out <- paste(cmd.out,"#SBATCH --time=168:00:00 --cpus-per-task=", num.cores,
    " --mem=64000M  --job-name=",paste(soup.sample.names[i], "presom", sep = "."),
    " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
    ### Add Read Group information
    cmd.out <- paste0(cmd.out, "java -jar ", picard, " AddOrReplaceReadGroups I=", file.path(soup.sample.dir.paths[i], bam.file),
        " RGID=1 RGPL=ILLUMINA RGPU=unit1 RGSM=20 RGLB=lib1 O=", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG.bam\n")))
    ### Mark Duplicate Reads
    cmd.out <- paste0(cmd.out, gatk, " MarkDuplicates -I ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG.bam")), " -O ",
        file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD.bam")), " -M ",
            file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD.txt")), "\n")
    ### Split NCigar Reads
    cmd.out <- paste0(cmd.out, gatk, " SplitNCigarReads -I ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD.bam")), " -O ",
        file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC.bam")), " -R ",
            opt$'ref-genome', "\n")
    ### Generate recalibration table based on covariates read group, reported quality score, machine cycle, and nucleotide context
    cmd.out <- paste0(cmd.out, gatk, " BaseRecalibrator -I ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC.bam")),
        " -R ", opt$'ref-genome', " --known-sites ", opt$'snp-db', " --known-sites ", opt$'indel-db', " -O ", 
            file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ1.table")), "\n")
    ### Run base quality score recalibration
    cmd.out <- paste0(cmd.out, gatk, " ApplyBQSR -R ", opt$'ref-genome', " -I ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC.bam")),
        " --bqsr-recal-file ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ1.table")),
            " -O ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ.bam")), "\n")
    ### Generate index file of BAM file
    cmd.out <- paste0(cmd.out, "java -jar ", picard, " BuildBamIndex I=", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ.bam")),
        "\n")
    ### Generate recalibration table on recalibrated base quality scores
    cmd.out <- paste0(cmd.out, gatk, " BaseRecalibrator -I ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ.bam")),
        " -R ", opt$'ref-genome', " --known-sites ", opt$'snp-db', " --known-sites ", opt$'indel-db', 
        " -O ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ2.table")), "\n")
    ### Analyze how recalibration worked across the covariates
    cmd.out <- paste0(cmd.out, gatk, " AnalyzeCovariates -before ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ1.table")),
        " -after ", file.path(sample.output.dir, paste0(soup.sample.names[i], "_RG_MD_SC_BQ2.table")), 
        " -plots ", file.path(sample.output.dir, paste0(soup.sample.names[i], ".pdf")), "\n")
    
    cat(cmd.out,file=script.filepath,append=F)

#### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
    if (!soup.sample.names[i] %in% completed.sample.names){
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }

## record what's so we can check for replicate sequencing data runs and append their reads together
    completed.sample.names <- append(completed.sample.names, soup.sample.names[i])       
}

system(paste("chmod 700 ",file.path(presom.script.dir,"presom.jobsub.bat"),sep=""))
#################
#completed.sample.names
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")




