#!/bin/bash

################ This script generates scripts to filter variants called by somatic variant caller GATK tool MuTect2 

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
    make_option(c("-v", "--vcf-dir"), type="character", default=NULL, 
              help="directory with subdirectories containing patient-level vcf files", metavar="character"),
    make_option(c("-b", "--biallelic-file"), type="character", default=NULL, 
              help="full directory to file with biallelic snps,", metavar="character"),
    make_option(c("-r", "--ref-genome"), type="character", default=NULL, 
              help="path to reference genome", metavar="character"),
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

filtersom.script.dir <- file.path(script.dir, "filtersom")
  if(file.exists(filtersom.script.dir)==F){
      cat("Generating somatic variant filtering script directory...\n")
      dir.create(filtersom.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(filtersom.script.dir,"filtersom.jobsub.bat")
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

    script.filepath <- file.path(filtersom.script.dir, paste(patient.id, "_som.sh", sep = ""))

################################### Print the script
    cmd.out <- NULL
    cmd.out <- paste("#!/bin/bash\n\n")
    cmd.out <- paste(cmd.out,"#SBATCH --time=168:00:00 --cpus-per-task=", num.cores,
    " --mem=64000M  --job-name=",paste(patient.id, "filtersom", sep = "."),
    " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
    cmd.out <- paste0(cmd.out, gatk, " GetPileupSummaries")
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

        ### continue with GetPileupSummaries options
        cmd.out <- paste0(cmd.out, " -V ", opt$'biallelic-file', " -L ", opt$'biallelic-file', " -O ", 
            file.path(patient.output.dir, paste0(patient.id, ".cond_pileups.table")), "\n")

        cmd.norm <- NULL

        ### check availability of normal samples
        if (cond.id  == "normal"){

            ### find samples corresponding to patient
            sample.ids <- anno.file$sample[anno.file$patient == patient.id & anno.file$condition == cond.id]

            ### loop through sample folders
            for (j in 1:length(sample.ids)){
                ### specify full path to sample
                sample.path <- soup.sample.dir.paths[grep(sample.ids[j], soup.sample.dir.paths)]

                ### specify full path to BAM file
                bam.file <- grep("_RG_MD_SC_BQ.bam", list.files(sample.path), value = TRUE)[1]

                ### add the string with full path to input BAM files for GetPileupSummaries
                cmd.norm <- paste0(cmd.norm, " -I ", file.path(sample.path, bam.file))

            }

            ### add GetPileupSummaries script for normal samples 
            cmd.out <- paste0(cmd.out, gatk, " GetPileupSummaries", cmd.norm, 
                " -V ", opt$'biallelic-file', " -L ", opt$'biallelic-file', " -O ",
                file.path(patient.output.dir, paste0(patient.id, ".normal_pileups.table")), "\n")
        }
    }

    ### add script for CalculateContamination to estimate contamination rate in each sample
    if (file.exists(file.path(patient.output.dir, ".normal_pileups.table"))){
        cmd.out <- paste0(cmd.out, gatk , " CalculateContamination -I", file.path(patient.output.dir, paste0(patient.id, ".cond_pileups.table")),
            " --matched-normal ", file.path(patient.output.dir, paste0(patient.id, ".normal_pileups.table")),
             " -O ", file.path(patient.output.dir, paste0(patient.id, ".contamination.table")), 
             " --tumor-segmentation ", file.path(patient.output.dir, paste0(patient.id, ".segments.table")), "\n") 
    } else {
        cmd.out <- paste0(cmd.out, gatk , " CalculateContamination -I", file.path(patient.output.dir, paste0(patient.id, ".cond_pileups.table")),
             " -O ", file.path(patient.output.dir, paste0(patient.id, ".contamination.table")), 
             " --tumor-segmentation ", file.path(patient.output.dir, paste0(patient.id, ".segments.table")), "\n\n")
    }

    ### specify full path to vcf file
    vcf.path <- file.path(opt$'vcf-dir', patient.id)
    vcf.file <- file.path(vcf.path, grep("vcf", list.files(file.path(opt$'vcf-dir', patient.id)), value = TRUE)[[1]][1])

    ### specify FilterMutectCalls
    cmd.out <- paste0(cmd.out, gatk , " FilterMutectCalls -R ", opt$'ref-genome', " -V ", vcf.file,
        " --contamination-table ", file.path(patient.output.dir, paste0(patient.id, ".contamination.table")),
        " --tumor-segmentation ", file.path(patient.output.dir, paste0(patient.id, ".segments.table")),
        " -O ", file.path(patient.output.dir, paste0(patient.id, "_filtered_mutect2.vcf\n\n")))


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
        
