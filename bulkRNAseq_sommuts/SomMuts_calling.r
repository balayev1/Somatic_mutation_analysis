#!/bin/bash

################ This script generates scripts to run somatic variant calling using GATK tool MuTect2 

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
    make_option(c("-s", "--germ-file"), type="character", default=NULL, 
              help="path to germline resource file ", metavar="character"),
    make_option(c("-d", "--pon"), type="character", default=NULL, 
              help="path to panel of normal file ", metavar="character"),
    make_option(c("-i", "--intervals"), type="character", default=NULL, 
              help="path to genome intervals file ", metavar="character"),
    make_option(c("-a", "--anno-file"), type="character", default=NULL, 
              help="path to tab-delimited sample annotation file", metavar="character")); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- opt$'parent-dir'
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}

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

som.script.dir <- file.path(script.dir, "som")
  if(file.exists(som.script.dir)==F){
      cat("Generating somatic variant caller script directory...\n")
      dir.create(som.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(som.script.dir,"som.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

### set the path to gatk executable
gatk <- "/home/abalay/scratch/Pipelines_and_Files/Packages/gatk/build/install/gatk/bin/gatk"

### set the path to picard executable
picard <- "/home/abalay/scratch/Pipelines_and_Files/Packages/picard/build/libs/picard-2.23.1-2-g3f22b39-SNAPSHOT-all.jar"

### load annotation table
anno.file <- read.table(opt$'anno-file', header=TRUE, sep="\t")
colnames(anno.file) <- c("sample", "patient", "condition")

## Extract the paths for the data parent directories (sample folders)
soup.sample.names <- list.files(input.parent.dir)
soup.sample.dir.paths <- file.path(input.parent.dir, soup.sample.names)

## make blank record for processed samples outside of the loop
completed.sample.names <- vector()

## Loop over patient folders
for(i in 1:length(unique(anno.file$patient))){

    patient.id <- unique(anno.file$patient)[i]

    patient.output.dir <- file.path(output.dir, patient.id)
    if (!dir.exists(patient.output.dir)){
        dir.create(patient.output.dir)
    }

    script.filepath <- file.path(som.script.dir, paste(patient.id, "_som.sh", sep = ""))

################################### Print the script
    cmd.out <- NULL
    cmd.out <- paste("#!/bin/bash\n\n")
    cmd.out <- paste(cmd.out,"#SBATCH --time=168:00:00 --cpus-per-task=", num.cores,
    " --mem=64000M  --job-name=",paste(patient.id, "som", sep = "."),
    " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
    cmd.out <- paste0(cmd.out, gatk, " Mutect2")
    ### verify number of conditions (max 2 allowed: one of conditions must be spelled as 'normal')
    for (k in 1:length(unique(anno.file$condition))){

        cond.id <- unique(anno.file$condition)[k]

        if (cond.id != "normal"){
            ### find samples corresponding to patient
            sample.ids <- anno.file$sample[anno.file$patient == patient.id & anno.file$condition == cond.id]
            
            ### loop through sample folders
            for (j in 1:length(sample.ids)){
                ### specify full path to sample
                sample.path <- soup.sample.dir.paths[grep(sample.ids[j], soup.sample.dir.paths)]

                ### specify full path to BAM file
                bam.file <- grep("_RG_MD_SC_BQ.bam", list.files(sample.path), value = TRUE)[1]


                ### add the string with full path to input BAM files for MuTect2
                cmd.out <- paste0(cmd.out, " -I ", file.path(sample.path, bam.file))
            }
        }

        cmd.norm <- NULL
        cmd.anno <- NULL

        if (cond.id  == "normal"){

            ### find samples corresponding to patient
            sample.ids <- anno.file$sample[anno.file$patient == patient.id & anno.file$condition == cond.id]

            ### loop through sample folders
            for (j in 1:length(sample.ids)){
                ### specify full path to sample
                sample.path <- soup.sample.dir.paths[grep(sample.ids[j], soup.sample.dir.paths)]

                ### specify full path to BAM file
                bam.file <- grep("_RG_MD_SC_BQ.bam", list.files(sample.path), value = TRUE)[1]

                ### add the string with full path to input BAM files for MuTect2
                cmd.norm <- paste0(cmd.norm, " -I ", file.path(sample.path, bam.file))

                ### add the string with sample annotation
                cmd.anno <- paste0(cmd.anno, " -normal ", strsplit(".bam", bam.file)[[1]][1])
            }
        }
    }

    ### continue with MuTect2 options
    cmd.out <- paste0(cmd.out, cmd.norm, cmd.anno)

    cmd.out <- paste0(cmd.out, " --germline-resource ", opt$'germ-file', " --panel-of-normals ", opt$'pon', 
        " --native-pair-hmm-threads ", num.cores, " --af-of-alleles-not-in-resource 0.0000025 ",
        "--reference ", opt$'ref-genome',
        " --genotype-germline-sites true ", 
        "--bam-output ", file.path(patient.output.dir, paste0(patient.id, "_unfiltered.bamout.bam ")),
        "--intervals ", opt$'intervals', " --interval-padding 100 ", 
        "-O ", file.path(patient.output.dir, paste0(patient.id, "_unfiltered_mutect2.vcf\n")))

    ### add MuTect2 script to the patient shell file
    cat(cmd.out,file=script.filepath,append=F)

    ### check if patient was processed and add it to job file
    if (!patient.id %in% completed.sample.names){
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)

    }

    completed.sample.names <- append(completed.sample.names, patient.id)
}

### make job file executable
system(paste0("chmod 700 ", jobsub.filepath))

#################
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")
        

