################# Filtering somatic mutations from RNA-Seq VCF files
################# Mutations called by MuTect2 software from GATK package

## request small number of resources
srun --mem=24GB --time=6:00:00 --pty --cpus-per-task=4 bash

## Activate conda environment
module load anaconda3
source activate r4_env

R

#### INFO on fields in VCF
##INFO=<ID=AS_FilterStatus,Number=A,Type=String,Description="Filter status for each allele, as assessed by ApplyRecalibration. Note that the VCF filter field will reflect the most lenient/sensitive sta$
##INFO=<ID=AS_SB_TABLE,Number=1,Type=String,Description="Allele-specific forward/reverse read counts for strand bias tests. Includes the reference and alleles separated by |.">
##INFO=<ID=AS_UNIQ_ALT_READ_COUNT,Number=A,Type=Integer,Description="Number of reads with unique start and mate end positions for each alt at a variant site">
##INFO=<ID=CONTQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to contamination">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
##INFO=<ID=ECNT,Number=1,Type=Integer,Description="Number of events in this haplotype">
##INFO=<ID=GERMQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not germline variants">
##INFO=<ID=MBQ,Number=R,Type=Integer,Description="median base quality">
##INFO=<ID=MFRL,Number=R,Type=Integer,Description="median fragment length">
##INFO=<ID=MMQ,Number=R,Type=Integer,Description="median mapping quality">
##INFO=<ID=MPOS,Number=A,Type=Integer,Description="median distance from end of read">
##INFO=<ID=NALOD,Number=A,Type=Float,Description="Negative log 10 odds of artifact in normal with same allele fraction as tumor">
##INFO=<ID=NCount,Number=1,Type=Integer,Description="Count of N bases in the pileup">
##INFO=<ID=NLOD,Number=A,Type=Float,Description="Normal log 10 likelihood ratio of diploid het or hom alt genotypes">
##INFO=<ID=OCM,Number=1,Type=Integer,Description="Number of alt reads whose original alignment doesn't match the current contig.">
##INFO=<ID=PON,Number=0,Type=Flag,Description="site found in panel of normals">
##INFO=<ID=POPAF,Number=A,Type=Float,Description="negative log 10 population allele frequencies of alt alleles">
##INFO=<ID=ROQ,Number=1,Type=Float,Description="Phred-scaled qualities that alt allele are not due to read orientation artifact">
##INFO=<ID=RPA,Number=R,Type=Integer,Description="Number of times tandem repeat unit is repeated, for each allele (including reference)">
##INFO=<ID=RU,Number=1,Type=String,Description="Tandem repeat unit (bases)">
##INFO=<ID=SEQQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles are not sequencing errors">
##INFO=<ID=STR,Number=0,Type=Flag,Description="Variant is a short tandem repeat">
##INFO=<ID=STRANDQ,Number=1,Type=Integer,Description="Phred-scaled quality of strand bias artifact">
##INFO=<ID=STRQ,Number=1,Type=Integer,Description="Phred-scaled quality that alt alleles in STRs are not polymerase slippage errors">
##INFO=<ID=TLOD,Number=A,Type=Float,Description="Log 10 likelihood ratio score of variant existing versus not existing">

## load required libraries
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)

## set the parent directory
input.dir <- "/home/abalay/scratch/Cell_Lines/SomMuts/SomMuts_filtering/sample_out"
if (!file.exists(input.dir)){
    cat("Input direcoty does not exist")
}

## set output directory
output.dir <- "/home/abalay/scratch/Cell_Lines/SomMuts"
if (!file.exists(output.dir)){
    cat("Generating output directory...")
    dir.create(output.dir)
}

figures.folder <- file.path(output.dir, "figures")
if (!file.exists(figures.folder)){
    cat("Generating figures folder...")
    dir.create(figures.folder)
}

setwd(figures.folder)

## list directories with vcf files
sample.names <- list.files(input.dir)
sample.dirs <- file.path(input.dir, sample.names)

## load VCF file and record number of each filterized variant category in each sample
germline <- c(); normal_artifact <- c(); panel_of_normals <- c(); strand_bias <- c(); slippage <- c(); orientation <- c()
contamination <- c(); pass <- c(); total_mutations <- c(); duplicate <- c(); low_allele_frac <- c()
base_qual <- c(); map_qual <- c(); multiallelic <- c(); fragment <- c(); position <- c(); 
clustered_events <- c(); weak_evidence <- c()

for (i in 1:length(sample.dirs)){
    ## set full path to vcf file
    vcf.path <- file.path(sample.dirs[i], grep("_filtered_mutect2.vcf", list.files(sample.dirs[i]), value=TRUE)[[1]][1])

    ## load the vcf file as tab-delimited
    vcf.file <- read.table(vcf.path, sep="\t", check.names=FALSE)

    ## define column names in vcf file
    colnames(vcf.file) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "20")

    total <- nrow(vcf.file) ## total # of called variants
    total_mutations <- append(total_mutations, total)

    ## fraction of variants found in germline resource
    germline <- append(germline, round(nrow(vcf.file[grep("germline", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found in matched controls
    normal_artifact <- append(normal_artifact, round(nrow(vcf.file[grep("normal_artifact", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found in panel of normals
    panel_of_normals <- append(panel_of_normals, round(nrow(vcf.file[grep("panel_of_normals", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to strand bias
    strand_bias <- append(strand_bias, round(nrow(vcf.file[grep("strand_bias", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to polymerase slippage
    slippage <- append(slippage, round(nrow(vcf.file[grep("slippage", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to library prep error
    orientation <- append(orientation, round(nrow(vcf.file[grep("orientation", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of contaminated variants
    contamination <- append(contamination, round(nrow(vcf.file[grep("contamination", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to duplicate insert start/end pairs of alt reads 
    duplicate <- append(duplicate, round(nrow(vcf.file[grep("duplicate", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to low allele fraction
    low_allele_frac <- append(low_allele_frac, round(nrow(vcf.file[grep("low_allele_frac", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to low base quality (median BQ of alt reads < 20)
    base_qual <- append(base_qual, round(nrow(vcf.file[grep("base_qual", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to low mapping quality (median MQ of alt reads < 30)
    map_qual <- append(map_qual, round(nrow(vcf.file[grep("map_qual", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to multiallelicity
    multiallelic <- append(multiallelic, round(nrow(vcf.file[grep("multiallelic", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to large difference in median fragment length between alt and ref reads
    fragment <- append(fragment, round(nrow(vcf.file[grep("fragment", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to position on the read
    position <- append(position, round(nrow(vcf.file[grep("position", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to clustered mutation on same haplotype (max < 2 mutations allowed)
    clustered_events <- append(clustered_events, round(nrow(vcf.file[grep("clustered_events", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants found due to potential sequencing error possibility
    weak_evidence <- append(weak_evidence, round(nrow(vcf.file[grep("weak_evidence", vcf.file$FILTER),])/total, 3) * 100)
    ## fraction of variants passing all filters
    pass <- append(pass, round(nrow(vcf.file[grep("PASS", vcf.file$FILTER),])/total, 3) * 100)
    cat("Sample ", sample.names[i], " is done\n")
}

## combine all info on categories to dataframe

df.combined <- data.frame(germline=germline, normal_artifact=normal_artifact, panel_of_normals=panel_of_normals,
    strand_bias=strand_bias, slippage=slippage, orientation=orientation, contamination=contamination, duplicate=duplicate,
    low_allele_frac=low_allele_frac, base_qual=base_qual, map_qual=map_qual, PASS=pass, clustered_events=clustered_events,
    multiallelic=multiallelic, fragment=fragment, position=position, weak_evidence=weak_evidence,
    row.names=sample.names)

total.df <- data.frame(Total=total_mutations, row.names=sample.names)

############################
########## Stacked barplot for each variant category
## remove all columns with colSums = 0
df.combined <- df.combined[,colSums(df.combined) != 0]

## make dataframe
samples <- c(); category <- c(); perc <- c()
for (i in 1:length(rownames(df.combined))){
    samples <- append(samples, rep(rownames(df.combined)[i], ncol(df.combined)))
    category <- append(category, colnames(df.combined))
    perc <- append(perc, as.numeric(df.combined[i,]))
}

df.stack <- data.frame(Sample = samples, Category = category, Percentage = perc)

## stacked bar plot
png("Stackedbarplot.huNSGmouse.celllines.053123.png", res=200, unit="in", height=8, width=14)
ggplot(df.stack, aes(fill=Category, y=Percentage, x=Sample)) + 
    geom_bar(position="stack", stat="identity")
dev.off()

## load RNA-editing DB
rna.editdb.path <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/TABLE1_hg38.txt"
rna.editdb.file <- fread(rna.editdb.path) 
rna.editdb.file <- paste0(rna.editdb.file$Region, rna.editdb.file$Position, rna.editdb.file$Ref, rna.editdb.file$Ed)

#########################################
######################## Keep only variants with filter status: PASS, clustered_events, fragment, position
for (j in 1:length(sample.dirs)){
    ## set full path to vcf file
    vcf.path <- file.path(sample.dirs[j], grep("_filtered_mutect2.vcf", list.files(sample.dirs[j]), value=TRUE)[[1]][1])

    ## load the vcf file as tab-delimited
    vcf.file <- read.table(vcf.path, sep="\t", check.names=FALSE)

    ## define column names in vcf file
    colnames(vcf.file) <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "20")

    ## keep only filters with PASS, fragment, position, clustered_events and haplotype
    filters.delete <- c("germline", "normal_artifact", "panel_of_normals", "strand_bias", "slippage", "contamination",
        "orientation", "duplicate", "low_allele_frac",
        "base_qual", "map_qual", "multiallelic", "weak_evidence")
    
    matches <- setdiff(seq(1, nrow(vcf.file)), grep(paste(filters.delete,collapse="|"), 
                        vcf.file$FILTER))
    
    vcf.file <- vcf.file[matches, ]

    ## subset variants with read depth >= 5
    variant.info <- c("1", "Fr|Rv", "DP", "ECNT", "GERMQ", "MBQ", "MFRL", "MMQ", "MPOS", "POPAF", "TLOD")
    temp.file <- vcf.file %>% separate(INFO, variant.info, ";")
    temp.file$DP <- as.numeric(gsub("DP=", "", temp.file$DP))
    dp.matches <- which(temp.file$DP >= 5)
    vcf.file <- vcf.file[dp.matches,]

    cat("Number of variants with DP >= 5 is", nrow(vcf.file), "\n")

    ## subset RNA-editing variants
    temp.file <- vcf.file
    temp.file <- paste0(temp.file$'#CHROM', temp.file$POS, temp.file$REF, temp.file$ALT)
    rna.edit.matches <- match(temp.file, rna.editdb.file)
    ## add rna editing sites to another file
    vcf.rnaedits <- vcf.file[which(!is.na(rna.edit.matches)), ]
    ## take non-rna editing sites
    vcf.file <- vcf.file[which(is.na(rna.edit.matches)), ]

    cat("Number of variants excluding RNA-editing variants = ", nrow(vcf.file), "\n")

    # save the main vcf file
    file=file.path(sample.dirs[j], paste0(sample.names[j], "_postfiltered_mutect2.vcf"))
    ## clear the file contents
    unlink(file)

    ### header
    vcf.lines <- readLines(vcf.path)
    write(vcf.lines[1:(grep("CHROM", vcf.lines))], file=file, append=TRUE)

    ## backbone
    write.table(vcf.file, file=file, append=TRUE, sep="\t", row.names=FALSE, quote=FALSE)

    ## save rna editing site vcf file
    rna.file=file.path(sample.dirs[j], paste0(sample.names[j], "_rnaedits.vcf"))
    ## clear the file contents
    unlink(rna.file)

    ### header
    write(vcf.lines[1:(grep("CHROM", vcf.lines))], file=rna.file, append=TRUE)

    ## backbone
    write.table(vcf.rnaedits, file=rna.file, append=TRUE, sep="\t", row.names=FALSE, quote=FALSE)

    cat("Clean VCF file for sample ", sample.names[j], " is done\n")
}










