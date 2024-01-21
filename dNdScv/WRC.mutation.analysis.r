####### Analysis of mutations at WRC/GYW sites


## Load required packages
library("ggpubr")
library("circlize")
library("ComplexHeatmap")
library("reshape2")


ebv.dndsout.path <- "/mnt/smb_share/AGS_AB/Cell_Lines/Mut_Profiles/dNdScv_IO/Output/EBV.huNSG.Sig_Mut.Genes@WRC.Result.txt"
ebv.kshv.dndsout.path <- "/mnt/smb_share/AGS_AB/Cell_Lines/Mut_Profiles/dNdScv_IO/Output/EBV+KSHV.huNSG.Sig_Mut.Genes@WRC.Result.txt"


## Load the gene list file with mutations
ebv.dndsout <- read.delim(ebv.dndsout.path, row.names=1, header=TRUE)
ebv.kshv.dndsout <- read.delim(ebv.kshv.dndsout.path, row.names=1, header=TRUE)

## Remove duplicated and unidentfied (i.e. no gene symbol) genes
gene.symbols <- gsub(".*:", "", rownames(ebv.dndsout)); ebv.dndsout <- ebv.dndsout[!duplicated(gene.symbols),]
gene.symbols <- gsub(".*:", "", rownames(ebv.kshv.dndsout)); ebv.kshv.dndsout <- ebv.kshv.dndsout[!duplicated(gene.symbols),]
rownames(ebv.dndsout) <- gsub(".*:", "", rownames(ebv.dndsout)); rownames(ebv.kshv.dndsout) <- gsub(".*:", "", rownames(ebv.kshv.dndsout))
ebv.dndsout <- ebv.dndsout[!rownames(ebv.dndsout) == "",]; ebv.kshv.dndsout <- ebv.kshv.dndsout[!rownames(ebv.kshv.dndsout) == "",]

colSums(ebv.dndsout)
# n_syn       n_mis       n_non       n_spl      n_utr5      n_utr3 
# 253.000     445.000     23.000     676.000     49.000     278.000
colSums(ebv.kshv.dndsout)
# n_syn       n_mis       n_non       n_spl      n_utr5      n_utr3 
# 216.000     335.000     29.000     768.000     43.000     215.000

#### Missense mutations
### Genes with pmis < 0.05
ebv.missense <- ebv.dndsout[ebv.dndsout$pmis_cv < 0.05, ]; ebv.missense <- ebv.missense[order(ebv.missense$pmis_cv),]
nrow(ebv.missense)
# [1] 205

ebv.kshv.missense <- ebv.kshv.dndsout[ebv.kshv.dndsout$pmis_cv < 0.05, ]; ebv.kshv.missense <- ebv.kshv.missense[order(ebv.kshv.missense$pmis_cv),]
nrow(ebv.kshv.missense)
# [1] 177

### Genes with >=2 missense mutations and missense rate > 1
ebv.mis.genes <- rownames(ebv.missense[ebv.missense$n_mis > 1 & ebv.missense$wmis_cv > 1,])
ebv.mis.genes
# [1] "ATP5MGL"  "RAB6C"    "RNF5"     "CCL4"     "IGLV1-47" "RPS17"   
#  [7] "UTP14C"   "TMPO"     "PABPC3"   "SLCO4C1"  "TUBB8B"   "IGHV3-30"
# [13] "AGAP5"    "IL4"      "GSTM5"    "ZSCAN22"  "PGLYRP1"  "NPR2"    
# [19] "MRPL45"   "SIGLEC11" "EPS8"     "DNAJC13"  "BRCA1"    "TAF1L"   
# [25] "IGHV1-69" "BCLAF3"   "TEP1"     "EHD1"     "CASD1"    "CELF6"   
# [31] "OCSTAMP"  "ACO2"     "TBC1D32"  "MYRF"     "IGHV3-74" "AMOT"    
# [37] "IGHV1-46" "FHOD3"    "MYOM1"    "ATRX"     "SCRIB"    "PLXNA3"  
# [43] "TRIOBP"   "VPS13C"

ebv.kshv.mis.genes <- rownames(ebv.kshv.missense[ebv.kshv.missense$n_mis > 1 & ebv.kshv.missense$wmis_cv > 1,])
ebv.kshv.mis.genes
# [1] "ATP5MGL"      "RAB6C"        "ERO1A"        "PTCD1"        "PABPC3"      
#  [6] "IGHV4-34"     "UTP14C"       "GSTM5"        "ATP5MF-PTCD1" "RPS17"       
# [11] "TMEM14C"      "IL4"          "HNRNPC"       "STEAP1B"      "TMPO"        
# [16] "LGMN"         "IGHV1-46"     "FAM114A1"     "EPS8"         "CELF6"       
# [21] "DTHD1"        "SOWAHB"       "FHIP2B"       "ATAD3B"       "DDX11"       
# [26] "RASA1"        "KIF4B"        "STAG3"        "RTEL1"        "CHD5"        
# [31] "DOCK7"        "ZNF831"       "ITPR3" 

intersect(ebv.mis.genes, ebv.kshv.mis.genes)
# [1] "ATP5MGL"  "RAB6C"    "RPS17"    "UTP14C"   "TMPO"     "PABPC3"  
#  [7] "IL4"      "GSTM5"    "EPS8"     "CELF6"    "IGHV1-46"

### Sample analysis of genes with missense mutations
## Load data
ebv.sample.dndsout <- read.delim("EBV.huNSG.Sample.Mut@WRC.Result.txt", header=TRUE); ebv.sample.dndsout$gene <- gsub(".*:", "", ebv.sample.dndsout$gene)
ebv.kshv.sample.dndsout <- read.delim("EBV+KSHV.huNSG.Sample.Mut@WRC.Result.txt", header=TRUE); ebv.kshv.sample.dndsout$gene <- gsub(".*:", "", ebv.kshv.sample.dndsout$gene)

## Annotate by infection status
ebv.sample.dndsout$infection <- rep("EBV", nrow(ebv.sample.dndsout)); ebv.kshv.sample.dndsout$infection <- rep("EBV+KSHV", nrow(ebv.kshv.sample.dndsout)); 

## Merge
sample.dndsout <- rbind(ebv.sample.dndsout, ebv.kshv.sample.dndsout)

######## Subset the genes with significant number of missense mutations
missense.genes <- c(setdiff(ebv.mis.genes, ebv.kshv.mis.genes), setdiff(ebv.kshv.mis.genes, ebv.mis.genes), intersect(ebv.mis.genes, ebv.kshv.mis.genes))
missense.genes
# [1] "RNF5"         "CCL4"         "IGLV1-47"     "SLCO4C1"      "TUBB8B"      
#  [6] "IGHV3-30"     "AGAP5"        "ZSCAN22"      "PGLYRP1"      "NPR2"        
# [11] "MRPL45"       "SIGLEC11"     "DNAJC13"      "BRCA1"        "TAF1L"       
# [16] "IGHV1-69"     "BCLAF3"       "TEP1"         "EHD1"         "CASD1"       
# [21] "OCSTAMP"      "ACO2"         "TBC1D32"      "MYRF"         "IGHV3-74"    
# [26] "AMOT"         "FHOD3"        "MYOM1"        "ATRX"         "SCRIB"       
# [31] "PLXNA3"       "TRIOBP"       "VPS13C"       "IGHV4-34"     "ERO1A"       
# [36] "PTCD1"        "ATP5MF-PTCD1" "TMEM14C"      "HNRNPC"       "STEAP1B"     
# [41] "LGMN"         "FAM114A1"     "DTHD1"        "FHIP2B"       "ATAD3B"      
# [46] "DDX11"        "SOWAHB"       "KIF4B"        "RASA1"        "STAG3"       
# [51] "RTEL1"        "CHD5"         "ZNF831"       "DOCK7"        "ITPR3"       
# [56] "ATP5MGL"      "RAB6C"        "RPS17"        "UTP14C"       "TMPO"        
# [61] "PABPC3"       "IL4"          "GSTM5"        "EPS8"         "CELF6"       
# [66] "IGHV1-46"

##### Dataframe: genes as rownames, # of coding, splice site, 5'UTR and 3'UTR mutations, # of samples in EBV and EBV+KSHV group
mutation.list <- list()
for (index in 1:length(missense.genes)){
    # Subset frame with gene of interest
    temp.dndsout <- sample.dndsout[sample.dndsout$gene == missense.genes[index], ]
    # Estimate number of mutations at coding, splice site, 5'UTR and 3'UTR of the gene
    n_missense <- nrow(temp.dndsout[temp.dndsout$impact == "Missense",])
    n_nonsense <- nrow(temp.dndsout[temp.dndsout$impact == "Nonsense",])
    n_splice <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_Splice",])
    n_5utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_5UTR",])
    n_3utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_3UTR",])
    # Estimate number of samples with mutated gene in EBV and EBV+KSHV group
    n_samp_ebv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV"]))
    n_sampl_ebv.kshv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV+KSHV"]))
    # Add all calculations into list
    temp.vector <- c(missense.genes[index],n_missense,n_nonsense,n_splice,n_5utr,n_3utr,n_samp_ebv,n_sampl_ebv.kshv)
    mutation.list[[index]] <- temp.vector
}

### Merge all information for genes
samp.mutation.object <- data.frame(do.call(rbind, mutation.list))

### Set column names for mutation object
colnames(samp.mutation.object) <- c("gene", "missense", "nonsense", "splice", "5utr", "3utr", "EBV-Lymphoma", "EBV+KSHV-Lymphoma")

### Check all missense genes are in object
all(missense.genes %in% samp.mutation.object$gene)
# [1] TRUE

### Set row names for mutation object
rownames(samp.mutation.object) <- samp.mutation.object[,1]; samp.mutation.object = samp.mutation.object[,!colnames(samp.mutation.object) %in% "gene"]

### Convert all columns from "character" to "numeric"
samp.mutation.object <- data.frame(apply(samp.mutation.object,2,as.numeric), check.names=FALSE, row.names=rownames(samp.mutation.object))

#################
##############
#### Make heatmap of mutations in genes per infection group
### Set colors
colors <- colorRamp2(c(0,1), c('white', 'blue'))

### Subset columns to make a heatmap
columns.keep <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")
mat <- as.matrix(samp.mutation.object[, colnames(samp.mutation.object) %in% columns.keep])

n_ebv_samp <- 6; n_ebv.kshv_samp <- 8 # number of samples with EBV and EBV+KSHV infections

## Convert matrix values to proportions
mat.prop <- as.matrix(cbind(mat[,1]/n_ebv_samp, mat[,2]/n_ebv.kshv_samp)); colnames(mat.prop) <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")

hp <- Heatmap(mat.prop, name = "Proportion of samples", col = colors, row_title = "Genes",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 7), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), 
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            if (mat.prop[i,j] > 0.8){
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(col = "white", fontsize = 10))
            }
            else if (mat.prop[i,j] < 0.8)
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#### Add stacked barplot of mutation types in genes
### Subset columns to make stacked barplot
columns.to.keep <- c("missense", "nonsense", "splice", "5utr", "3utr")
mut.type.object <- samp.mutation.object[, colnames(samp.mutation.object) %in% columns.to.keep]

mut.type.object <- t(apply(mut.type.object, 1, function(x) x/sum(x)))

### Set colors for stacked barplot annotation
color.anno <- c("red", "blue", "green", "purple", "orange")
sp <- rowAnnotation(Fraction = anno_barplot(mut.type.object, gp = gpar(fill = color.anno), bar_width = 1, width = unit(4, "cm"),
    axis_param = list(at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1),side="bottom"), axis = TRUE))


#### Add legend for stacked barplot
lgd = Legend(labels = colnames(mut.type.object), title = "Mutations", legend_gp = gpar(col = color.anno), 
    direction="vertical", type = "points", pch = 16)


png("Heatmap+stacked_barplot.genes.missense.mutations.WRC.png", res=200, unit="in", height=8, width=11)
draw(hp + sp, heatmap_legend_side="bottom", annotation_legend_list = lgd, annotation_legend_side="right")
dev.off()


#### Splice site and Nonsense mutations
#### Nonsense mutations:
### Genes with ptrunc < 0.05
ebv.trunc <- ebv.dndsout[ebv.dndsout$ptrunc_cv < 0.05, ]; ebv.trunc <- ebv.trunc[order(ebv.trunc$ptrunc_cv),]
nrow(ebv.trunc)
# [1] 457

ebv.kshv.trunc <- ebv.kshv.dndsout[ebv.kshv.dndsout$ptrunc_cv < 0.05, ]; ebv.kshv.trunc <- ebv.kshv.trunc[order(ebv.kshv.trunc$ptrunc_cv),]
nrow(ebv.kshv.trunc)
# [1] 504


### Genes with at least 1 nonsense mutation and nonsense rate > 1
ebv.non.genes <- rownames(ebv.trunc[ebv.trunc$n_non > 0 & ebv.trunc$wnon_cv > 1,])
ebv.non.genes
#  [1] "BIRC3"        "NUP205"       "PIBF1"        "FAM102A"      "PDIA5"       
#  [6] "ATF6B"        "DACT3"        "SENP3-EIF4A1" "CCDC86"       "OXA1L"       
# [11] "ABCD3"        "PEX6"         "AJM1"         "MCM9"         "ANKRD27"     
# [16] "NUTM1"        "SBNO2"        "MYO1B"        "DAAM2"        "TIMELESS"    
# [21] "SMC1B"

ebv.kshv.non.genes <- rownames(ebv.kshv.trunc[ebv.kshv.trunc$n_non > 0 & ebv.kshv.trunc$wnon_cv > 1,])
ebv.kshv.non.genes
# [1] "VPS50"      "MYRF"       "PTGS1"      "BCORL1"     "NRXN3"     
#  [6] "SECISBP2L"  "ST6GALNAC4" "POLA1"      "MYO18B"     "RAB12"     
# [11] "NAXD"       "SLC25A22"   "ANKRD9"     "CNOT6"      "SPIRE2"    
# [16] "SLC25A25"   "WDR11"      "LMF1"       "PLCL1"      "PPFIA3"    
# [21] "TUBGCP3"    "ERCC6"      "CDC42BPG"   "CRYBG3"     "DNAJC13"   
# [26] "TSC22D1"    "SZT2"

intersect(ebv.non.genes, ebv.kshv.non.genes)
# character(0)

######## Subset the genes with significant number of missense mutations
nonsense.genes <- c(ebv.non.genes, ebv.kshv.non.genes)
nonsense.genes
# [1] "BIRC3"        "NUP205"       "PIBF1"        "FAM102A"      "PDIA5"       
#  [6] "ATF6B"        "DACT3"        "SENP3-EIF4A1" "CCDC86"       "OXA1L"       
# [11] "ABCD3"        "PEX6"         "AJM1"         "MCM9"         "ANKRD27"     
# [16] "NUTM1"        "SBNO2"        "MYO1B"        "DAAM2"        "TIMELESS"    
# [21] "SMC1B"        "VPS50"        "MYRF"         "PTGS1"        "BCORL1"      
# [26] "NRXN3"        "SECISBP2L"    "ST6GALNAC4"   "POLA1"        "MYO18B"      
# [31] "RAB12"        "NAXD"         "SLC25A22"     "ANKRD9"       "CNOT6"       
# [36] "SPIRE2"       "SLC25A25"     "WDR11"        "LMF1"         "PLCL1"       
# [41] "PPFIA3"       "TUBGCP3"      "ERCC6"        "CDC42BPG"     "CRYBG3"      
# [46] "DNAJC13"      "TSC22D1"      "SZT2"   

##### Dataframe: genes as rownames, # of coding, splice site, 5'UTR and 3'UTR mutations, # of samples in EBV and EBV+KSHV group
mutation.list <- list()
for (index in 1:length(nonsense.genes)){
    # Subset frame with gene of interest
    temp.dndsout <- sample.dndsout[sample.dndsout$gene == nonsense.genes[index], ]
    # Estimate number of mutations at coding, splice site, 5'UTR and 3'UTR of the gene
    n_missense <- nrow(temp.dndsout[temp.dndsout$impact == "Missense",])
    n_nonsense <- nrow(temp.dndsout[temp.dndsout$impact == "Nonsense",])
    n_splice <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_Splice",])
    n_5utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_5UTR",])
    n_3utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_3UTR",])
    # Estimate number of samples with mutated gene in EBV and EBV+KSHV group
    n_samp_ebv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV"]))
    n_sampl_ebv.kshv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV+KSHV"]))
    # Add all calculations into list
    temp.vector <- c(nonsense.genes[index],n_missense,n_nonsense,n_splice,n_5utr,n_3utr,n_samp_ebv,n_sampl_ebv.kshv)
    mutation.list[[index]] <- temp.vector
}

### Merge all information for genes
samp.mutation.object <- data.frame(do.call(rbind, mutation.list))

### Set column names for mutation object
colnames(samp.mutation.object) <- c("gene", "missense", "nonsense", "splice", "5utr", "3utr", "EBV-Lymphoma", "EBV+KSHV-Lymphoma")

### Check all missense genes are in object
all(nonsense.genes %in% samp.mutation.object$gene)
# [1] TRUE

### Set row names for mutation object
rownames(samp.mutation.object) <- samp.mutation.object[,1]; samp.mutation.object = samp.mutation.object[,!colnames(samp.mutation.object) %in% "gene"]

### Convert all columns from "character" to "numeric"
samp.mutation.object <- data.frame(apply(samp.mutation.object,2,as.numeric), check.names=FALSE, row.names=rownames(samp.mutation.object))

#################
##############
#### Make heatmap of mutations in genes per infection group
### Set colors
colors <- colorRamp2(c(0,1), c('white', 'blue'))

### Subset columns to make a heatmap
columns.keep <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")
mat <- as.matrix(samp.mutation.object[, colnames(samp.mutation.object) %in% columns.keep])

n_ebv_samp <- 6; n_ebv.kshv_samp <- 8 # number of samples with EBV and EBV+KSHV infections

## Convert matrix values to proportions
mat.prop <- as.matrix(cbind(mat[,1]/n_ebv_samp, mat[,2]/n_ebv.kshv_samp)); colnames(mat.prop) <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")

hp <- Heatmap(mat.prop, name = "Proportion of samples", col = colors, row_title = "Genes",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 7), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), 
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            if (mat.prop[i,j] > 0.8){
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(col = "white", fontsize = 10))
            }
            else if (mat.prop[i,j] < 0.8)
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#### Add stacked barplot of mutation types in genes
### Subset columns to make stacked barplot
columns.to.keep <- c("missense", "nonsense", "splice", "5utr", "3utr")
mut.type.object <- samp.mutation.object[, colnames(samp.mutation.object) %in% columns.to.keep]

mut.type.object <- t(apply(mut.type.object, 1, function(x) x/sum(x)))

### Set colors for stacked barplot annotation
color.anno <- c("red", "blue", "green", "purple", "orange")
sp <- rowAnnotation(Fraction = anno_barplot(mut.type.object, gp = gpar(fill = color.anno), bar_width = 1, width = unit(4, "cm"),
    axis_param = list(at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1),side="bottom"), axis = TRUE))


#### Add legend for stacked barplot
lgd = Legend(labels = colnames(mut.type.object), title = "Mutations", legend_gp = gpar(col = color.anno), 
    direction="vertical", type = "points", pch = 16)


png("Heatmap+stacked_barplot.genes.nonsense.mutations.WRC.png", res=200, unit="in", height=8, width=11)
draw(hp + sp, heatmap_legend_side="bottom", annotation_legend_list = lgd, annotation_legend_side="right")
dev.off()





#### Splice site mutations:
### Genes with > 3 splice site mutations and splice site rate > 1
ebv.spl.genes <- rownames(ebv.trunc[ebv.trunc$n_spl > 3 & ebv.trunc$wspl_cv > 1,])
ebv.spl.genes
# [1] "SLC35A5"    "RNF14"      "FOCAD"      "NIN"        "MARK3"     
#  [6] "ZNF141"     "CCNJ"       "VPS13D"     "ST7"        "FNBP4"     
# [11] "EP300"      "URB1"       "CCP110"     "KMT2C"      "ITGA1"     
# [16] "KIF20A"     "POMT1"      "CEP131"     "LRCH3"      "INTS4"     
# [21] "GLDC"       "MYO18B"     "KIF15"      "LY75"       "LY75-CD302"
# [26] "VPS13A"     "MCM10"

ebv.kshv.spl.genes <- rownames(ebv.kshv.trunc[ebv.kshv.trunc$n_spl > 3 & ebv.kshv.trunc$wspl_cv > 1,])
ebv.kshv.spl.genes
# [1] "SLC35B1" "SLC35A5" "CCNJ"    "MARK3"   "CEP131"  "NFAT5"   "SRD5A3" 
#  [8] "ZC3H7B"  "CNTRL"   "OTUD6B"  "PI4KA"   "EP300"   "URB1"    "VPS13A" 
# [15] "FURIN"   "SPATS2"  "VPS13D"  "NIN"     "HCFC2"   "SAFB"    "RBM27"  
# [22] "FNBP4"   "SLU7"    "USP47"   "KMT2C"   "FOCAD"

intersect(ebv.spl.genes, ebv.kshv.spl.genes)
# [1] "SLC35A5" "FOCAD"   "NIN"     "MARK3"   "CCNJ"    "VPS13D"  "FNBP4"  
#  [8] "EP300"   "URB1"    "KMT2C"   "CEP131"  "VPS13A" 

######## Subset the genes with significant number of splice site mutations
splice.genes <- c(setdiff(ebv.spl.genes, ebv.kshv.spl.genes), setdiff(ebv.kshv.spl.genes,ebv.spl.genes),intersect(ebv.spl.genes, ebv.kshv.spl.genes))
splice.genes
# [1] "RNF14"      "ZNF141"     "ST7"        "CCP110"     "ITGA1"     
#  [6] "KIF20A"     "POMT1"      "LRCH3"      "INTS4"      "GLDC"      
# [11] "MYO18B"     "KIF15"      "LY75"       "LY75-CD302" "MCM10"     
# [16] "SLC35B1"    "NFAT5"      "SRD5A3"     "ZC3H7B"     "CNTRL"     
# [21] "OTUD6B"     "PI4KA"      "FURIN"      "SPATS2"     "HCFC2"     
# [26] "SAFB"       "RBM27"      "SLU7"       "USP47"      "SLC35A5"   
# [31] "FOCAD"      "NIN"        "MARK3"      "CCNJ"       "VPS13D"    
# [36] "FNBP4"      "EP300"      "URB1"       "KMT2C"      "CEP131"    
# [41] "VPS13A" 

##### Dataframe: genes as rownames, # of coding, splice site, 5'UTR and 3'UTR mutations, # of samples in EBV and EBV+KSHV group
mutation.list <- list()
for (index in 1:length(splice.genes)){
    # Subset frame with gene of interest
    temp.dndsout <- sample.dndsout[sample.dndsout$gene == splice.genes[index], ]
    # Estimate number of mutations at coding, splice site, 5'UTR and 3'UTR of the gene
    n_missense <- nrow(temp.dndsout[temp.dndsout$impact == "Missense",])
    n_nonsense <- nrow(temp.dndsout[temp.dndsout$impact == "Nonsense",])
    n_splice <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_Splice",])
    n_5utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_5UTR",])
    n_3utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_3UTR",])
    # Estimate number of samples with mutated gene in EBV and EBV+KSHV group
    n_samp_ebv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV"]))
    n_sampl_ebv.kshv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV+KSHV"]))
    # Add all calculations into list
    temp.vector <- c(splice.genes[index],n_missense,n_nonsense,n_splice,n_5utr,n_3utr,n_samp_ebv,n_sampl_ebv.kshv)
    mutation.list[[index]] <- temp.vector
}

### Merge all information for genes
samp.mutation.object <- data.frame(do.call(rbind, mutation.list))

### Set column names for mutation object
colnames(samp.mutation.object) <- c("gene", "missense", "nonsense", "splice", "5utr", "3utr", "EBV-Lymphoma", "EBV+KSHV-Lymphoma")

### Check all missense genes are in object
all(splice.genes %in% samp.mutation.object$gene)
# [1] TRUE

### Set row names for mutation object
rownames(samp.mutation.object) <- samp.mutation.object[,1]; samp.mutation.object = samp.mutation.object[,!colnames(samp.mutation.object) %in% "gene"]

### Convert all columns from "character" to "numeric"
samp.mutation.object <- data.frame(apply(samp.mutation.object,2,as.numeric), check.names=FALSE, row.names=rownames(samp.mutation.object))

#################
##############
#### Make heatmap of mutations in genes per infection group
### Set colors
colors <- colorRamp2(c(0,1), c('white', 'blue'))

### Subset columns to make a heatmap
columns.keep <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")
mat <- as.matrix(samp.mutation.object[, colnames(samp.mutation.object) %in% columns.keep])

n_ebv_samp <- 6; n_ebv.kshv_samp <- 8 # number of samples with EBV and EBV+KSHV infections

## Convert matrix values to proportions
mat.prop <- as.matrix(cbind(mat[,1]/n_ebv_samp, mat[,2]/n_ebv.kshv_samp)); colnames(mat.prop) <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")

hp <- Heatmap(mat.prop, name = "Proportion of samples", col = colors, row_title = "Genes",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 7), 
    column_names_gp = gpar(fontsize = 8),
    heatmap_legend_param = list(direction="horizontal"), 
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            if (mat.prop[i,j] > 0.8){
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(col = "white", fontsize = 10))
            }
            else if (mat.prop[i,j] < 0.8)
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#### Add stacked barplot of mutation types in genes
### Subset columns to make stacked barplot
columns.to.keep <- c("missense", "nonsense", "splice", "5utr", "3utr")
mut.type.object <- samp.mutation.object[, colnames(samp.mutation.object) %in% columns.to.keep]

mut.type.object <- t(apply(mut.type.object, 1, function(x) x/sum(x)))

### Set colors for stacked barplot annotation
color.anno <- c("red", "blue", "green", "purple", "orange")
sp <- rowAnnotation(Fraction = anno_barplot(mut.type.object, gp = gpar(fill = color.anno), bar_width = 1, width = unit(4, "cm"),
    axis_param = list(at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1),side="bottom"), axis = TRUE))


#### Add legend for stacked barplot
lgd = Legend(labels = colnames(mut.type.object), title = "Mutations", legend_gp = gpar(col = color.anno), 
    direction="vertical", type = "points", pch = 16)


png("Heatmap+stacked_barplot.genes.splice.site.mutations.WRC.png", res=200, unit="in", height=8, width=11)
draw(hp + sp, heatmap_legend_side="bottom", annotation_legend_list = lgd, annotation_legend_side="right")
dev.off()




#### 5'UTR mutations
### Genes with > 0 5'UTR mutations and 5'UTR rate > 1
ebv.utr5 <- ebv.dndsout[ebv.dndsout$n_utr5 > 0 & ebv.dndsout$wutr5_cv > 1, ]
nrow(ebv.utr5)
# [1] 22

ebv.utr5.genes <- rownames(ebv.utr5)
ebv.utr5.genes
# [1] "MECP2"     "MTRNR2L8"  "FHL2"      "MTRNR2L12" "CXCL2"     "ZNF334"   
#  [7] "PSMA2"     "TAGAP"     "HPSE"      "GAB2"      "CPT1B"     "SOGA1"    
# [13] "GTSF1L"    "CLIC1"     "RC3H2"     "ACD"       "MTRNR2L4"  "ZBED2"    
# [19] "ZNF232"    "CAMK2D"    "ZGLP1"     "RCOR2"

ebv.kshv.utr5 <- ebv.kshv.dndsout[ebv.kshv.dndsout$n_utr5 > 0 & ebv.kshv.dndsout$wutr5_cv > 1, ]
nrow(ebv.kshv.utr5)
# [1] 20

ebv.kshv.utr5.genes <- rownames(ebv.kshv.utr5)
ebv.kshv.utr5.genes
# [1] "MTRNR2L8"      "MTRNR2L12"     "USP51"         "TOX"          
#  [5] "IL12RB1"       "SYNJ2BP"       "SYNJ2BP-COX16" "TVP23C-CDRT4" 
#  [9] "ASB9"          "ZNF432"        "ZNF425"        "FRK"          
# [13] "GMEB2"         "SOS2"          "PHKA1"         "DCLRE1A"      
# [17] "LRRC25"        "RPA3"          "KLF9"          "LEF1"  

intersect(ebv.utr5.genes, ebv.kshv.utr5.genes)
# [1] "MTRNR2L8"  "MTRNR2L12"

######## Subset the genes with significant number of 5'UTR mutations
utr5.genes <- c(setdiff(ebv.utr5.genes, ebv.kshv.utr5.genes), setdiff(ebv.kshv.utr5.genes,ebv.utr5.genes),intersect(ebv.utr5.genes, ebv.kshv.utr5.genes))
utr5.genes
# [1] "MECP2"         "FHL2"          "CXCL2"         "ZNF334"       
#  [5] "PSMA2"         "TAGAP"         "HPSE"          "GAB2"         
#  [9] "CPT1B"         "SOGA1"         "GTSF1L"        "CLIC1"        
# [13] "RC3H2"         "ACD"           "MTRNR2L4"      "ZBED2"        
# [17] "ZNF232"        "CAMK2D"        "ZGLP1"         "RCOR2"        
# [21] "USP51"         "TOX"           "IL12RB1"       "SYNJ2BP"      
# [25] "SYNJ2BP-COX16" "TVP23C-CDRT4"  "ASB9"          "ZNF432"       
# [29] "ZNF425"        "FRK"           "GMEB2"         "SOS2"         
# [33] "PHKA1"         "DCLRE1A"       "LRRC25"        "RPA3"         
# [37] "KLF9"          "LEF1"          "MTRNR2L8"      "MTRNR2L12"

##### Dataframe: genes as rownames, # of coding, splice site, 5'UTR and 3'UTR mutations, # of samples in EBV and EBV+KSHV group
mutation.list <- list()
for (index in 1:length(utr5.genes)){
    # Subset frame with gene of interest
    temp.dndsout <- sample.dndsout[sample.dndsout$gene == utr5.genes[index], ]
    # Estimate number of mutations at coding, splice site, 5'UTR and 3'UTR of the gene
    n_missense <- nrow(temp.dndsout[temp.dndsout$impact == "Missense",])
    n_nonsense <- nrow(temp.dndsout[temp.dndsout$impact == "Nonsense",])
    n_splice <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_Splice",])
    n_5utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_5UTR",])
    n_3utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_3UTR",])
    # Estimate number of samples with mutated gene in EBV and EBV+KSHV group
    n_samp_ebv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV"]))
    n_sampl_ebv.kshv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV+KSHV"]))
    # Add all calculations into list
    temp.vector <- c(utr5.genes[index],n_missense,n_nonsense,n_splice,n_5utr,n_3utr,n_samp_ebv,n_sampl_ebv.kshv)
    mutation.list[[index]] <- temp.vector
}

### Merge all information for genes
samp.mutation.object <- data.frame(do.call(rbind, mutation.list))

### Set column names for mutation object
colnames(samp.mutation.object) <- c("gene", "missense", "nonsense", "splice", "5utr", "3utr", "EBV-Lymphoma", "EBV+KSHV-Lymphoma")

### Check all missense genes are in object
all(utr5.genes %in% samp.mutation.object$gene)
# [1] TRUE

### Set row names for mutation object
rownames(samp.mutation.object) <- samp.mutation.object[,1]; samp.mutation.object = samp.mutation.object[,!colnames(samp.mutation.object) %in% "gene"]

### Convert all columns from "character" to "numeric"
samp.mutation.object <- data.frame(apply(samp.mutation.object,2,as.numeric), check.names=FALSE, row.names=rownames(samp.mutation.object))

#################
##############
#### Make heatmap of mutations in genes per infection group
### Set colors
colors <- colorRamp2(c(0,1), c('white', 'blue'))

### Subset columns to make a heatmap
columns.keep <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")
mat <- as.matrix(samp.mutation.object[, colnames(samp.mutation.object) %in% columns.keep])

n_ebv_samp <- 6; n_ebv.kshv_samp <- 8 # number of samples with EBV and EBV+KSHV infections

## Convert matrix values to proportions
mat.prop <- as.matrix(cbind(mat[,1]/n_ebv_samp, mat[,2]/n_ebv.kshv_samp)); colnames(mat.prop) <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")

hp <- Heatmap(mat.prop, name = "Proportion of samples", col = colors, row_title = "Genes",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 7), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), 
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            if (mat.prop[i,j] > 0.8){
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(col = "white", fontsize = 10))
            }
            else if (mat.prop[i,j] < 0.8)
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#### Add stacked barplot of mutation types in genes
### Subset columns to make stacked barplot
columns.to.keep <- c("missense", "nonsense", "splice", "5utr", "3utr")
mut.type.object <- samp.mutation.object[, colnames(samp.mutation.object) %in% columns.to.keep]

mut.type.object <- t(apply(mut.type.object, 1, function(x) x/sum(x)))

### Set colors for stacked barplot annotation
color.anno <- c("red", "blue", "green", "purple", "orange")
sp <- rowAnnotation(Fraction = anno_barplot(mut.type.object, gp = gpar(fill = color.anno), bar_width = 1, width = unit(4, "cm"),
    axis_param = list(at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1),side="bottom"), axis = TRUE))


#### Add legend for stacked barplot
lgd = Legend(labels = colnames(mut.type.object), title = "Mutations", legend_gp = gpar(col = color.anno), 
    direction="vertical", type = "points", pch = 16)


png("Heatmap+stacked_barplot.genes.5UTR.mutations.WRC.png", res=200, unit="in", height=8, width=11)
draw(hp + sp, heatmap_legend_side="bottom", annotation_legend_list = lgd, annotation_legend_side="right")
dev.off()


#### 3'UTR mutations
### Genes with > 1 3'UTR mutations and 3'UTR rate > 1
ebv.utr3 <- ebv.dndsout[ebv.dndsout$n_utr3 > 1 & ebv.dndsout$wutr3_cv > 1, ]
nrow(ebv.utr3)
# [1] 35

ebv.utr3.genes <- rownames(ebv.utr3)
ebv.utr3.genes
# [1] "UQCRHL"   "H3-5"     "MACC1"    "C10orf95" "FAM3A"    "ZNF607"  
#  [7] "ARHGAP4"  "POTEC"    "BCLAF3"   "SOX5"     "AKT2"     "XRCC2"   
# [13] "MBOAT1"   "ZNF558"   "ESR2"     "CFAP418"  "ZNF529"   "TTC33"   
# [19] "DESI1"    "ZNHIT6"   "DENND6B"  "SVIP"     "NPTXR"    "BRWD3"   
# [25] "SLC25A53" "SETD7"    "CLOCK"    "SLC7A11"  "DAPK2"    "ADAM10"  
# [31] "RPL37"    "CSNK2A1"  "CCDC127"  "STX7"     "INO80D"

ebv.kshv.utr3 <- ebv.kshv.dndsout[ebv.kshv.dndsout$n_utr3 > 1 & ebv.kshv.dndsout$wutr3_cv > 1, ]
nrow(ebv.kshv.utr3)
# [1] 22

ebv.kshv.utr3.genes <- rownames(ebv.kshv.utr3)
ebv.kshv.utr3.genes
# [1] "UQCRHL"   "H3-5"     "ZNF607"   "FAM3A"    "C10orf95" "CTBP2"   
#  [7] "SLN"      "LOXL2"    "ATRX"     "C1QTNF6"  "KLHL8"    "NEK4"    
# [13] "AVPR1A"   "RAB14"    "SORT1"    "SVIP"     "MACC1"    "SYN3"    
# [19] "IKZF3"    "NEGR1"    "CSNK2A1"  "SOD2"  

intersect(ebv.utr3.genes, ebv.kshv.utr3.genes)
# [1] "UQCRHL"   "H3-5"     "MACC1"    "C10orf95" "FAM3A"    "ZNF607"   "SVIP"    
# [8] "CSNK2A1"

######## Subset the genes with significant number of 3'UTR mutations
utr3.genes <- c(setdiff(ebv.utr3.genes, ebv.kshv.utr3.genes), setdiff(ebv.kshv.utr3.genes,ebv.utr3.genes),intersect(ebv.utr3.genes, ebv.kshv.utr3.genes))
utr3.genes
# [1] "ARHGAP4"  "POTEC"    "BCLAF3"   "SOX5"     "AKT2"     "XRCC2"   
#  [7] "MBOAT1"   "ZNF558"   "ESR2"     "CFAP418"  "ZNF529"   "TTC33"   
# [13] "DESI1"    "ZNHIT6"   "DENND6B"  "NPTXR"    "BRWD3"    "SLC25A53"
# [19] "SETD7"    "CLOCK"    "SLC7A11"  "DAPK2"    "ADAM10"   "RPL37"   
# [25] "CCDC127"  "STX7"     "INO80D"   "CTBP2"    "SLN"      "LOXL2"   
# [31] "ATRX"     "C1QTNF6"  "KLHL8"    "NEK4"     "AVPR1A"   "RAB14"   
# [37] "SORT1"    "SYN3"     "IKZF3"    "NEGR1"    "SOD2"     "UQCRHL"  
# [43] "H3-5"     "MACC1"    "C10orf95" "FAM3A"    "ZNF607"   "SVIP"    
# [49] "CSNK2A1" 


##### Dataframe: genes as rownames, # of coding, splice site, 5'UTR and 3'UTR mutations, # of samples in EBV and EBV+KSHV group
mutation.list <- list()
for (index in 1:length(utr3.genes)){
    # Subset frame with gene of interest
    temp.dndsout <- sample.dndsout[sample.dndsout$gene == utr3.genes[index], ]
    # Estimate number of mutations at coding, splice site, 5'UTR and 3'UTR of the gene
    n_missense <- nrow(temp.dndsout[temp.dndsout$impact == "Missense",])
    n_nonsense <- nrow(temp.dndsout[temp.dndsout$impact == "Nonsense",])
    n_splice <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_Splice",])
    n_5utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_5UTR",])
    n_3utr <- nrow(temp.dndsout[temp.dndsout$impact == "Essential_3UTR",])
    # Estimate number of samples with mutated gene in EBV and EBV+KSHV group
    n_samp_ebv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV"]))
    n_sampl_ebv.kshv <- length(unique(temp.dndsout$sampleID[temp.dndsout$infection == "EBV+KSHV"]))
    # Add all calculations into list
    temp.vector <- c(utr3.genes[index],n_missense,n_nonsense,n_splice,n_5utr,n_3utr,n_samp_ebv,n_sampl_ebv.kshv)
    mutation.list[[index]] <- temp.vector
}

### Merge all information for genes
samp.mutation.object <- data.frame(do.call(rbind, mutation.list))

### Set column names for mutation object
colnames(samp.mutation.object) <- c("gene", "missense", "nonsense", "splice", "5utr", "3utr", "EBV-Lymphoma", "EBV+KSHV-Lymphoma")

### Check all missense genes are in object
all(utr3.genes %in% samp.mutation.object$gene)
# [1] TRUE

### Set row names for mutation object
rownames(samp.mutation.object) <- samp.mutation.object[,1]; samp.mutation.object = samp.mutation.object[,!colnames(samp.mutation.object) %in% "gene"]

### Convert all columns from "character" to "numeric"
samp.mutation.object <- data.frame(apply(samp.mutation.object,2,as.numeric), check.names=FALSE, row.names=rownames(samp.mutation.object))

#################
##############
#### Make heatmap of mutations in genes per infection group
### Set colors
colors <- colorRamp2(c(0,1), c('white', 'blue'))

### Subset columns to make a heatmap
columns.keep <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")
mat <- as.matrix(samp.mutation.object[, colnames(samp.mutation.object) %in% columns.keep])

n_ebv_samp <- 6; n_ebv.kshv_samp <- 8 # number of samples with EBV and EBV+KSHV infections

## Convert matrix values to proportions
mat.prop <- as.matrix(cbind(mat[,1]/n_ebv_samp, mat[,2]/n_ebv.kshv_samp)); colnames(mat.prop) <- c("EBV-Lymphoma", "EBV+KSHV-Lymphoma")

hp <- Heatmap(mat.prop, name = "Proportion of samples", col = colors, row_title = "Genes",
    row_title_rot = 90, column_names_rot = 45, cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = "left", 
    column_names_side = "top", row_names_gp = gpar(fontsize = 7), 
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(direction="horizontal"), 
    cell_fun = function(j, i, x, y, width, height,fill) 
        {
            if (mat.prop[i,j] > 0.8){
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(col = "white", fontsize = 10))
            }
            else if (mat.prop[i,j] < 0.8)
                grid::grid.text(sprintf("%.f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
        })

#### Add stacked barplot of mutation types in genes
### Subset columns to make stacked barplot
columns.to.keep <- c("missense", "nonsense", "splice", "5utr", "3utr")
mut.type.object <- samp.mutation.object[, colnames(samp.mutation.object) %in% columns.to.keep]

mut.type.object <- t(apply(mut.type.object, 1, function(x) x/sum(x)))

### Set colors for stacked barplot annotation
color.anno <- c("red", "blue", "green", "purple", "orange")
sp <- rowAnnotation(Fraction = anno_barplot(mut.type.object, gp = gpar(fill = color.anno), bar_width = 1, width = unit(4, "cm"),
    axis_param = list(at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1),side="bottom"), axis = TRUE))


#### Add legend for stacked barplot
lgd = Legend(labels = colnames(mut.type.object), title = "Mutations", legend_gp = gpar(col = color.anno), 
    direction="vertical", type = "points", pch = 16)


png("Heatmap+stacked_barplot.genes.3UTR.mutations.WRC.png", res=200, unit="in", height=8, width=11)
draw(hp + sp, heatmap_legend_side="bottom", annotation_legend_list = lgd, annotation_legend_side="right")
dev.off()


#####################
##################
#### Make stacked barplot: x-axis: type of mutations by genomic location, y-axis: number of genes with certain number of \
#### mutations at respective genomic location
### Define mutation types
types <- colnames(ebv.dndsout)[2:6]
types
# [1] "n_mis"  "n_non"  "n_spl"  "n_utr5" "n_utr3"

### Full names of mutation types
full.names <- c("Missense", "Nonsense", "Splice site", "5'UTR", "3'UTR")

####### Dataframe: rownames: mutation type + number, ebv: number of genes in EBV group, ebv+kshv: number of genes in EBV+KSHV group
temp.ebv.vector <- c(); temp.ebv.kshv.vector <- c(); temp.anno.vector <- c(); mutation.type.vector <- c()
for (type in 1:length(types)) {
    sel_cv <- paste("w", substr(types[type],3,nchar(types[type])), "_cv", sep="")
    quantity <- max(max(ebv.dndsout[,types[type]]), max(ebv.kshv.dndsout[,types[type]]))
    # Quantify number of genes per mutation type in EBV and EBV+KSHV group of samples
    for (i in 1:quantity) {
        n.ebv <- nrow(ebv.dndsout[ebv.dndsout[types[type]] == i & ebv.dndsout[sel_cv] > 1,])
        temp.ebv.vector <- append(temp.ebv.vector, n.ebv)

        n.ebv.kshv <- nrow(ebv.kshv.dndsout[ebv.kshv.dndsout[types[type]] == i & ebv.kshv.dndsout[sel_cv] > 1,])
        temp.ebv.kshv.vector <- append(temp.ebv.kshv.vector, n.ebv.kshv)

        temp.anno.vector <- append(temp.anno.vector, paste(substr(sel_cv,2,5), "_", i, sep=""))

        # Add mutation type column
        mutation.type.vector <- append(mutation.type.vector, full.names[type])
    }
}

## Dataframe
num.genes.type <- rbind(temp.ebv.vector, temp.ebv.kshv.vector)
colnames(num.genes.type) <- temp.anno.vector; rownames(num.genes.type) <- c("EBV", "EBV+KSHV")
num.genes.type <- data.frame(t(num.genes.type))

## Add mutation type column
num.genes.type$mut.type <- mutation.type.vector

## Remove columns with all zeros
num.genes.type <- num.genes.type[rowSums(num.genes.type[,-ncol(num.genes.type)]) != 0, ]

## Add number of mutations column
num.genes.type$mut.number <- as.character(lapply(rownames(num.genes.type), function(x) gsub(".*_","", x)))
num.genes.type$mut.number <- factor(num.genes.type$mut.number, levels=1:16)

### Grouped stacked barplot
melted <- melt(num.genes.type, c("mut.type", "mut.number"))

melted$group <- ''
melted[melted$variable == 'EBV',]$group <- "EBV"
melted[melted$variable != 'EBV',]$group <- "EBV+KSHV"

png("Grouped_stacked_barplot.numgenes.permuttype.WRC.png", res=200, unit="in", height=8, width=11)
ggplot(melted, aes(x = group, y = value, fill = mut.number, label = value)) + 
    geom_bar(stat = 'identity', position = 'stack') + 
    facet_grid(~ mut.type) +
    theme(legend.position = "right") +
    scale_fill_discrete(name = "Number of mutations") +
    ylab("Number of genes") +
    ggtitle("Number of genes per mutation type by number")
dev.off()


####### Genes as left rows, Infection as top column, metric: raw number of mutations and proportion, stacked barplot on right row
####### 
gene.count.file <- "/mnt/smb_share/AGS_AB/Cell_Lines/Transcriptomic_Profiles/Count.Matrices/huNSGTissue_gene_count.csv"
gene.count.matrix <- read.delim(gene.count.file, header=TRUE)


















