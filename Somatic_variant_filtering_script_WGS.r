################# Filtering somatic mutations from Whole Genome Sequencing (WGS) VCF files 
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

parser <- ArgumentParser(description='Process some integers')
parser$add_argument('--work-dir', dest='workdir', type="character", nargs=1, help='Working directory to vcf files')
parser$add_argument('--anno-dir', dest='annodir', type="character", nargs=1, help='Directory to annotation file')

args <- parser$parse_args()

dirs <- list.dirs(args$workdir)
dirs <- dirs[grepl(paste(c("05-", "06-", "07-", "09-", "HTMCP-"),collapse="|"), dirs)]
dirs
# [1] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/05-12939"         
#  [2] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/05-25439"         
#  [3] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/05-25674"         
#  [4] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-11535"         
#  [5] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-15256"         
#  [6] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-16716"         
#  [7] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-19919"         
#  [8] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-22057"         
#  [9] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/06-34043"         
# [10] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/07-35482"         
# [11] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/09-33003"         
# [12] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/HTMCP-01-01-00003"
# [13] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/HTMCP-01-01-00012"
# [14] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/HTMCP-01-02-00013"
# [15] "/mnt/smb_share/AGS_AB/WGS_RAW_DLBCL/SNVFiles/HTMCP-01-02-00017"


### Upload VCF samples and record number of mutations in each filtering category + PASS
germline <- c(); normal_artifact <- c(); panel_of_normals <- c(); strand_bias <- c(); slippage <- c(); orientation <- c()
contamination <- c(); pass <- c(); sample_name <- c(); total_mutations <- c(); duplicate <- c(); low_allele_frac <- c()
strict_strand <- c(); base_qual <- c(); map_qual <- c()
for (i in 1:length(dirs)){
    ### Load and convert VCF file into dataframe
    vcf.list <- list.files(dirs[i], pattern = "_filtered_mutect2.vcf", full.names = TRUE)
    vcf.file <- read.table(vcf.list[1], sep="\t", check.names=FALSE, quote="")

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
    strand_bias=strand_bias, slippage=slippage, orientation=orientation, contamination=contamination, PASS=pass,
    row.names=sample_name)

total.df <- data.frame(Total=total_mutations, row.names=sample_name)

#################
##############
#### Make heatmap of % mutations for each filter
### Set colors
colors <- colorRamp2(c(0,1), c('blue', 'red'))

### Make barplot
bp <- HeatmapAnnotation(Number_of_variants = anno_barplot(total.df$Total, which="column", axis = TRUE, 
    gp = gpar(fill="green", fontsize = 4), bar_width = 1))

hp <- Heatmap(t(as.matrix(df.combined)), name = "Proportion of variants", col = colors, row_title = "Filtering Category",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 10), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), top_annotation = bp,
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            grid::grid.text(sprintf("%.3f", t(as.matrix(df.combined))[i,j]), x, y, gp = gpar(col = "white", fontsize = 10))
        })

png("Heatmap.filtering.percentage.WGS.DLBCL.png", res=200, unit="in", height=8, width=14)
draw(hp, heatmap_legend_side="right")
dev.off()


###########
###### Remove specific filtered variants
anno.file <- read_excel(paste(args$annodir,"/EBV-transformed_lymphoblast_cell_lines.xlsx", sep=""), sheet = "DBGAP WGS DLBCL NORMAL")  # excel annotation file for WGS data

for (i in 1:length(dirs)){
    cat("Starting processing sample ", basename(dirs[i]), "\n")
    ### Load and convert VCF file into dataframe
    vcf.list <- list.files(dirs[i], pattern = "_filtered_mutect2.vcf", full.names = TRUE)
    vcf.file <- read.table(vcf.list[1], sep="\t", check.names=FALSE)

    vcf.lines <- readLines(vcf.list[1])
    # vcf.file <- vcf.file[grep("CHROM", vcf.lines):length(vcf.lines),]
    colnames(vcf.file) <- strsplit(vcf.lines[grep("CHROM", vcf.lines)], "\t")[[1]]

    if (length(which(is.na(vcf.file$FILTER))) > 0){
        vcf.file <- vcf.file[-(which(is.na(vcf.file$FILTER))),]
    }

    ######### Remove all variants with below specified filters
    filter.remove <- c("panel_of_normals", "strand_bias", "slippage", "orientation", "contamination", 
    "duplicate", "low_allele_frac", "strict_strand")
    for (j in 1:length(filter.remove)) {
        if (length(grep(filter.remove[j], vcf.file$FILTER)) != 0){
            vcf.file <- vcf.file[-(grep(filter.remove[j], vcf.file$FILTER)),] 
        }
    }

    ######### Variants that are tumor in normal contamination
    ### Subset all germline and normal_artifact files
    germline.indices <- grep("germline", vcf.file$FILTER); normal_artifact.indices <- grep("normal_artifact", vcf.file$FILTER)
    indices <- c(setdiff(germline.indices, normal_artifact.indices), intersect(germline.indices, normal_artifact.indices),
        setdiff(normal_artifact.indices, germline.indices))
    vcf.germ <- vcf.file[indices,]

    ### Tumor in normal variants
    ### Denote tumor and normal samples
    tumor <- c(); normal <- c()
    for (index in 1:length(anno.file$`SRA-ID`)){
        if (anno.file$`Tumor (yes/no)`[index] == "yes"){
            tumor <- append(tumor, paste(anno.file$`SRA-ID`[index], "_dbGaP-28434", sep=""))
        }
        if (anno.file$`Tumor (yes/no)`[index] == "no"){
            normal <- append(normal, paste(anno.file$`SRA-ID`[index], "_dbGaP-28434", sep=""))
        }
    }

    ### Extract allelic depths of alternative alleles in tumor samples
    alt.allele.count.tumor <- rep(0, nrow(vcf.germ))
    total.allele.count.tumor <- rep(0, nrow(vcf.germ))
    for (k in tumor){
        if (k %in% colnames(vcf.germ)) {
            temp <- sapply(strsplit(vcf.germ[,k], split=":"), "[[", 2)
            ad.alt <- as.numeric(as.vector(sapply(strsplit(temp, split=","), "[[", 2)))
            ad.ref <- as.numeric(as.vector(sapply(strsplit(temp, split=","), "[", 1)))
            ad.total <- ad.alt + ad.ref
            alt.allele.count.tumor <- alt.allele.count.tumor + ad.alt
            total.allele.count.tumor <- total.allele.count.tumor + ad.total
        }
    }

    af.tumor <- sapply(alt.allele.count.tumor/total.allele.count.tumor, function(x) if (!is.na(x)) round(x, 3) else x=0)

    ### Extract allelic depths of alternative alleles in normal samples
    alt.allele.count.normal <- rep(0, nrow(vcf.germ))
    total.allele.count.normal <- rep(0, nrow(vcf.germ))
    for (k in normal){
        if (k %in% colnames(vcf.germ)) {
            temp <- sapply(strsplit(vcf.germ[,k], split=":"), "[[", 2)
            ad.alt <- as.numeric(as.vector(sapply(strsplit(temp, split=","), "[[", 2)))
            ad.ref <- as.numeric(as.vector(sapply(strsplit(temp, split=","), "[", 1)))
            ad.total <- ad.alt + ad.ref
            alt.allele.count.normal <- alt.allele.count.normal + ad.alt
            total.allele.count.normal <- total.allele.count.normal + ad.total
        }
    }

    af.normal <- sapply(alt.allele.count.normal/total.allele.count.normal, function(x) if (!is.na(x)) round(x, 3) else x=0)

    ### Dataframe of allele frequency in tumor and normal samples
    germ.df <- data.frame(AF_Tumor = af.tumor, AF_Normal = af.normal)


    ### Plot AF in tumor vs normal
    png(paste(basename(dirs[i]), ".Tumor_vs_Normal.germline.png", sep=""), res=200, unit="in", height=8, width=14)
    print({
        p <- ggplot(germ.df, aes(x=log10(AF_Normal), y=AF_Tumor)) + ggtitle("Allele frequency of variants \n labeled as germline") +
        xlab("AF in Normal (log10)") + ylab("AF in DLBCL") + geom_point() + theme_minimal() + 
        geom_vline(xintercept = log10(0.0005), colour="red", linetype = "longdash") + 
        geom_hline(yintercept = 0.5, colour="red", linetype = "longdash")
        p})
    dev.off()


    ### Remove all variants not passing germline and normal_artifact filters with AFs in tumor < 0.5 and in normal > log10(0.0005)
    kept.indices <- which(germ.df$AF_Tumor > 0.5 & log10(germ.df$AF_Normal) < log10(0.0005))

    # % of kept variants (tumor in normal contamination)
    cat("Fraction of kept variants as tumor in normal contamination = ", length(kept.indices)/nrow(germ.df), "\n")

    vcf.germ <- vcf.germ[kept.indices,]
    vcf.file <- rbind(vcf.germ, vcf.file[-(indices),])



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
    system(paste("sudo cp ", file, " ", dirs[i], sep = ""))
    cat("Sample ", basename(dirs[i]), " is done\n")
}








