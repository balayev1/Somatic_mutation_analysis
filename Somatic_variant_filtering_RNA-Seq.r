################# Filtering somatic mutations from RNA-Seq VCF tumor-only files for cell lines
################# Mutations called by MuTect2 software from GATK package

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

### Load library
library("ComplexHeatmap")
library("circlize")
library("readxl")
library("ggplot2")
library("argparse")

######## Step 1: Heatmap of % mutations for each sample labeled in following filtered categories: germline, normal_artifact, 
######## panel_of_normals, strand_bias, slippage, orientation, contamination + PASS 

parser <- ArgumentParser(description='Somatic variant filtering for cell lines')
parser$add_argument('--work-dir', dest='workdir', type="character", nargs=1, help='Working directory to vcf files')

args <- parser$parse_args()

dirs <- list.dirs(args$workdir)
dirs <- dirs[grepl(paste(c("-R1", "-R2", "-R3"),collapse="|"), dirs)]
dirs
# [1] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/AP2-R1"   
#  [2] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/AP3-R1"   
#  [3] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/AP5-R1"   
#  [4] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/BCBL1-R1" 
#  [5] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E1-R1"    
#  [6] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E2-R1"    
#  [7] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E3-R1"    
#  [8] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-3-R1"  
#  [9] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-3-R2"  
# [10] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-3-R3"  
# [11] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-9-R1"  
# [12] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-9-R2"  
# [13] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-9-R3"  
# [14] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E4-R1"    
# [15] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E5-R1"    
# [16] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E6-R1"    
# [17] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/E7-R1"    
# [18] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK3-11-R1"
# [19] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK3-11-R2"
# [20] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK3-11-R3"
# [21] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK3-13-R1"
# [22] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK3-13-R2"
# [23] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK4-11-R1"
# [24] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK4-11-R2"
# [25] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/EK4-11-R3"
# [26] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/HBL6-R1"  
# [27] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL1-R1"  
# [28] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL1-R2"  
# [29] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL1-R3"  
# [30] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL2-R1"  
# [31] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL2-R2"  
# [32] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL2-R3"  
# [33] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL3-R1"  
# [34] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL3-R2"  
# [35] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/LCL3-R3"  
# [36] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/M1-R1"    
# [37] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/M2-R1"    
# [38] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/M3-R1"    
# [39] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/M4-R1"    
# [40] "/home/cluster/abalay/scratch/Cell_Lines/SNVFiles/M5-R1"


### Upload VCF samples and record number of mutations in each filtering category + PASS
germline <- c(); normal_artifact <- c(); panel_of_normals <- c(); strand_bias <- c(); slippage <- c(); orientation <- c()
contamination <- c(); pass <- c(); sample_name <- c(); total_mutations <- c(); duplicate <- c(); low_allele_frac <- c()
strict_strand <- c(); base_qual <- c(); map_qual <- c()
for (i in 1:length(dirs)){
    ### Load and convert VCF file into dataframe
    vcf.list <- list.files(dirs[i], pattern = "_filtered_mutect2.vcf", full.names = TRUE)
    vcf.file <- read.table(vcf.list[1], sep="\t", check.names=FALSE)

    vcf.lines <- readLines(vcf.list[1])
    vcf.file <- vcf.file[grep("CHROM", vcf.lines):length(vcf.lines),]
    colnames(vcf.file) <- strsplit(vcf.lines[grep("CHROM", vcf.lines)], "\t")[[1]]
    ### Estimate % of each category 
    sample_name <- append(sample_name, basename(dirs[i]))
    total <- nrow(vcf.file)
    total_mutations <- append(total_mutations, total)
    germline <- append(germline, round(nrow(vcf.file[grep("germline", vcf.file$FILTER),])/total, 3))
    normal_artifact <- append(normal_artifact, round(nrow(vcf.file[grep("normal_artifact", vcf.file$FILTER),])/total, 3))
    panel_of_normals <- append(panel_of_normals, round(nrow(vcf.file[grep("panel_of_normals", vcf.file$FILTER),])/total, 3))
    strand_bias <- append(strand_bias, round(nrow(vcf.file[grep("strand_bias", vcf.file$FILTER),])/total, 3))
    slippage <- append(slippage, round(nrow(vcf.file[grep("slippage", vcf.file$FILTER),])/total, 3))
    orientation <- append(orientation, round(nrow(vcf.file[grep("orientation", vcf.file$FILTER),])/total, 3))
    contamination <- append(contamination, round(nrow(vcf.file[grep("contamination", vcf.file$FILTER),])/total, 3))
    duplicate <- append(duplicate, round(nrow(vcf.file[grep("duplicate", vcf.file$FILTER),])/total, 3))
    low_allele_frac <- append(low_allele_frac, round(nrow(vcf.file[grep("low_allele_frac", vcf.file$FILTER),])/total, 3))
    strict_strand <- append(strict_strand, round(nrow(vcf.file[grep("strict_strand", vcf.file$FILTER),])/total, 3))
    base_qual <- append(base_qual, round(nrow(vcf.file[grep("base_qual", vcf.file$FILTER),])/total, 3))
    map_qual <- append(map_qual, round(nrow(vcf.file[grep("map_qual", vcf.file$FILTER),])/total, 3))
    pass <- append(pass, round(nrow(vcf.file[grep("PASS", vcf.file$FILTER),])/total, 3))
    cat("Sample ", basename(dirs[i]), " is done\n")
}

## Combine all info on categories to dataframe
df.combined <- data.frame(germline=germline, normal_artifact=normal_artifact, panel_of_normals=panel_of_normals,
    strand_bias=strand_bias, slippage=slippage, orientation=orientation, contamination=contamination, duplicate=duplicate,
    low_allele_frac=low_allele_frac, strict_strand=strict_strand, base_qual=base_qual, map_qual=map_qual, PASS=pass,
    row.names=sample_name)

total.df <- data.frame(Total=total_mutations, row.names=sample_name)

#################
##############
#### Make heatmap of % mutations for each filter
#### Heatmap 1: huNSG mice cell lines
cell.lines1 <- c("AP2-R1", "AP3-R1", "AP5-R1", "BCBL1-R1", "HBL6-R1", "E4-3-R1", "E4-3-R2", "E4-3-R3", "E4-9-R1", "E4-9-R2",
    "E4-9-R3", "EK3-11-R1", "EK3-11-R2", "EK3-11-R3", "EK3-13-R1", "EK3-13-R2", "EK4-11-R1", "EK4-11-R2", "EK4-11-R3",
    "LCL1-R1", "LCL1-R2", "LCL1-R3", "LCL2-R1", "LCL2-R2", "LCL2-R3", "LCL3-R1", "LCL3-R2", "LCL3-R3")
df1.combined <- df.combined[rownames(df.combined) %in% cell.lines1,]
total.df$Sample <- rownames(total.df); total.df1 <-  total.df[total.df$Sample %in% cell.lines1,]

### Set colors
colors <- colorRamp2(c(0,1), c('blue', 'red'))

### Make barplot
bp <- HeatmapAnnotation(Number_of_variants = anno_barplot(total.df1$Total, which="column", axis = TRUE, 
    gp = gpar(fill="green", fontsize = 4), bar_width = 1))

hp <- Heatmap(t(as.matrix(df1.combined)), name = "Proportion of variants", col = colors, row_title = "Filtering Category",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 10), rect_gp = gpar(col = "red", lwd = 1), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), top_annotation = bp,
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            grid::grid.text(sprintf("%.3f", t(as.matrix(df1.combined))[i,j]), x, y, gp = gpar(col = "white", fontsize = 10))
        })

png("Heatmap.filtering.percentage.Cell.LinesWT.png", res=200, unit="in", height=8, width=14)
draw(hp, heatmap_legend_side="right")
dev.off()


#### Heatmap 2: huNSG mice cell lines HLA2-
cell.lines2 <- c("E1-R1", "E2-R1", "E3-R1", "E4-R1", "E5-R1", "E6-R1", "E7-R1", "M1-R1", "M2-R1", "M3-R1", "M4-R1", "M5-R1")
df2.combined <- df.combined[rownames(df.combined) %in% cell.lines2,]
total.df2 <-  total.df[total.df$Sample %in% cell.lines2,]

### Set colors
colors <- colorRamp2(c(0,1), c('blue', 'red'))

### Make barplot
bp <- HeatmapAnnotation(Number_of_variants = anno_barplot(total.df2$Total, which="column", axis = TRUE, 
    gp = gpar(fill="green", fontsize = 4), bar_width = 1))

hp <- Heatmap(t(as.matrix(df2.combined)), name = "Proportion of variants", col = colors, row_title = "Filtering Category",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 10), rect_gp = gpar(col = "red", lwd = 1), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), top_annotation = bp,
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            grid::grid.text(sprintf("%.3f", t(as.matrix(df2.combined))[i,j]), x, y, gp = gpar(col = "white", fontsize = 10))
        })

png("Heatmap.filtering.percentage.Cell.LinesHLA2-.png", res=200, unit="in", height=8, width=14)
draw(hp, heatmap_legend_side="right")
dev.off()

###########
###### Remove specific filtered variants
for (i in 1:length(dirs)){
    cat("Starting processing sample ", basename(dirs[i]), "\n")
    ### Load and convert VCF file into dataframe
    vcf.list <- list.files(dirs[i], pattern = "_filtered_mutect2.vcf", full.names = TRUE)
    vcf.file <- read.table(vcf.list[1], sep="\t", check.names=FALSE, quote="")

    vcf.lines <- readLines(vcf.list[1])
    colnames(vcf.file) <- strsplit(vcf.lines[grep("CHROM", vcf.lines)], "\t")[[1]]

    if (length(which(is.na(vcf.file$FILTER))) > 0){
        vcf.file <- vcf.file[-(which(is.na(vcf.file$FILTER))),]
    }

    ######### Remove all variants with below specified filters
    filter.remove <- c("panel_of_normals", "strand_bias", "slippage", "orientation", "contamination", "germline", "normal_artifact", 
        "duplicate", "low_allele_frac", "strict_strand")
    for (j in 1:length(filter.remove)) {
        if (length(grep(filter.remove[j], vcf.file$FILTER)) != 0){
            vcf.file <- vcf.file[-(grep(filter.remove[j], vcf.file$FILTER)),] 
        }
    }

    ############# Remove variants with low median base quality and median mapping quality
    ### Subset all above mentioned files
    mbq.values <- sapply(strsplit(vcf.file$INFO, split=";"), "[", 6)
    mbq.values <- as.numeric(as.vector(sapply(strsplit(mbq.values, split=","), "[", 2)))

    mmq.values <- sapply(strsplit(vcf.file$INFO, split=";"), "[", 8)
    mmq.values <- as.numeric(as.vector(sapply(strsplit(mmq.values, split=","), "[", 2)))

    ### Dataframe of median base qualities and mapping qualities
    qual.df <- data.frame(MBQ = mbq.values, MMQ = mmq.values)

    ### Plot MBQ vs MMQ
    png(paste(basename(dirs[i]), ".MBQvsMMQ.png", sep=""), res=200, unit="in", height=8, width=14)
    print({
        p <- ggplot(qual.df, aes(x=MMQ, y=MBQ)) + ggtitle("MBQ VS MMQ of all variants") +
        xlab("MMQ") + ylab("MBQ") + geom_point() + theme_minimal() + 
        geom_vline(xintercept = 20, colour="red", linetype = "longdash") + 
        geom_hline(yintercept = 15, colour="red", linetype = "longdash")
        p})
    dev.off()

    mbq.indices <- grep("base_qual", vcf.file$FILTER)
    mmq.indices <- grep("map_qual", vcf.file$FILTER)
    indices <- c(setdiff(mbq.indices, mmq.indices), intersect(mbq.indices, mmq.indices),
        setdiff(mmq.indices, mbq.indices))

    qual.df$status <- rep("PASS", nrow(qual.df))
    qual.df$status[indices] = "FAIL"

    ### Plot MBQ vs MMQ by PASS or FAIL
    png(paste(basename(dirs[i]), ".MBQvsMMQbyStatus.png", sep=""), res=200, unit="in", height=8, width=14)
    print({
        p <- ggplot(qual.df, aes(x=MMQ, y=MBQ, color=status)) + ggtitle("MBQ VS MMQ of all variants") +
        xlab("MMQ") + ylab("MBQ") + geom_point() + theme_minimal() + scale_color_manual(values=c("blue", "red")) +
        geom_vline(xintercept = 20, colour="red", linetype = "longdash") + 
        geom_hline(yintercept = 15, colour="red", linetype = "longdash")
        p})
    dev.off()

    kept.indices <- which(qual.df$MBQ > 15 & qual.df$MMQ > 20)

    # % of kept variants after base and mapping quality observation
    cat("Fraction of kept variants after base_qual and map_qual filtering = ", length(kept.indices)/nrow(qual.df), "\n")

    vcf.file <- vcf.file[kept.indices,]


    ################# Filter any variant labeled as weak_evidence
    if (length(grep("weak_evidence", vcf.file$FILTER)) != 0){
            vcf.file <- vcf.file[-(grep("weak_evidence", vcf.file$FILTER)),] 
    }
    
    cat("Number of post-filtering variants = ", nrow(vcf.file), "\n")
    
    ### Order by chromosomes
    vcf.file[,"#CHROM"]<-factor(vcf.file[,"#CHROM"], levels=unique(vcf.file[,"#CHROM"]))
    vcf.file <- vcf.file[order(vcf.file[,"#CHROM"], vcf.file[,"POS"]),]

    ### Save to file
    file=paste(basename(dirs[i]), "_post_filtered_mutect2.vcf", sep ="")
    write(vcf.lines[1:(grep("CHROM", vcf.lines)-1)], file=file, append=TRUE)
    write.table(vcf.file, file=file, append=TRUE, sep="\t", row.names=FALSE, quote=FALSE)
    system(paste("cp ", file, " ", dirs[i], sep = ""))
    cat("Sample ", basename(dirs[i]), " is done\n")
}


















