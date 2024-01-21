#' @param cdsfile Path to the reference transcript table.
#' @param genomefile Path to the indexed reference genome file.
#' @param outfile Output file name (default = "RefCDS.rda").
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate
#' @param excludechrs Vector or string with chromosome names to be excluded from the RefCDS object (default: no chromosome will be excluded). The mitochondrial chromosome should be excluded as it has different genetic code and mutation rates, either using the excludechrs argument or not including mitochondrial transcripts in cdsfile.
#' @param onlychrs Vector of valid chromosome names (default: all chromosomes will be included)
#' @param useids Combine gene IDs and gene names (columns 1 and 2 of the input table) as long gene names (default = F)
#' 
#' @export

buildref = function(cdsfile, genomefile, outfile = "RefCDS.rda", numcode = 1, excludechrs = NULL, onlychrs = NULL, useids = F) {
    
    ## 1. Valid chromosomes and reference CDS per gene
    message("[1/3] Preparing the environment...")
    
    reftable = read.table(cdsfile, header=1, sep="\t", stringsAsFactors=F, quote="\"", na.strings="-", fill=TRUE) # Loading the reference table
    colnames(reftable) = c("gene.id","gene.name","cds.id","chr","chr.coding.start","chr.coding.end","cds.start","cds.end","length","strand","tss","5utr.start", "5utr.end", "3utr.start", "3utr.end")
    reftable[,5:15] = suppressWarnings(lapply(reftable[,5:15], as.numeric)) # Convert to numeric
    
    # Checking for systematic absence of gene names (it happens in some BioMart inputs)
    longname = paste(reftable$gene.id, reftable$gene.name, sep=":") # Gene name combining the Gene stable ID and the Associated gene name
    if (useids==T) {
        reftable$gene.name = longname # Replacing gene names by a concatenation of gene ID and gene name
    }
    if (length(unique(reftable$gene.name))<length(unique(longname))) {
        warning(sprintf("%0.0f unique gene IDs (column 1) found. %0.0f unique gene names (column 2) found. Consider combining gene names and gene IDs or replacing gene names by gene IDs to avoid losing genes (see useids argument in ? buildref)",length(unique(reftable$gene.id)),length(unique(reftable$gene.name))))
    }
    
    # Reading chromosome names in the fasta file using its index file if it already exists or creating it if it does not exist. The index file is also used by scanFa later.
    validchrs = as.character(GenomicRanges::seqnames(Rsamtools::scanFaIndex(genomefile)))

    validchrs = setdiff(validchrs, excludechrs)
    if (length(onlychrs)>0) {
        validchrs = validchrs[validchrs %in% onlychrs]
    }
    
    # Restricting to chromosomes present in both the genome file and the CDS table
    if (any(validchrs %in% unique(reftable$chr))) {
        validchrs = validchrs[validchrs %in% unique(reftable$chr)]
    } else { # Try adding a chr prefix
        reftable$chr = paste("chr", reftable$chr, sep="")
        validchrs = validchrs[validchrs %in% unique(reftable$chr)]
        if (length(validchrs)==0) { # No matching chromosome names
            stop("No chromosome names in common between the genome file and the CDS table")
        }
    }

# Removing genes that fall partially or completely outside of the available chromosomes/contigs
    
    reftable = reftable[reftable[,1]!="" & reftable[,2]!="" & reftable[,3]!="" & !is.na(reftable[,5]) & !is.na(reftable[,6]),] # Removing invalid entries
    reftable = reftable[which(reftable$chr %in% validchrs),] # Only valid chromosomes
    
    transc_gr = GenomicRanges::GRanges(reftable$chr, IRanges::IRanges(reftable$chr.coding.start,reftable$chr.coding.end))
    chrs_gr = Rsamtools::scanFaIndex(genomefile)
    ol = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="within", select="all"))
    
    # Issuing an error if any transcript falls outside of the limits of a chromosome. Possibly due to a mismatch between the assembly used for the reference table and the reference genome.
    if (length(unique(ol[,1])) < nrow(reftable)) {
        stop(sprintf("Aborting buildref. %0.0f rows in cdsfile have coordinates that fall outside of the corresponding chromosome length. Please ensure that you are using the same assembly for the cdsfile and genomefile",nrow(reftable)-length(unique(ol[,1]))))
    }
    
    reftable = reftable[unique(ol[,1]),] 
    
    # Identifying genes starting or ending at the ends of a chromosome/contig
    # Because buildref and dndscv need to access the base before and after each coding position, genes overlapping the ends
    # of a contig will be trimmed by three bases and a warning will be issued listing those genes.
    
    fullcds = intersect(reftable$cds.id[reftable$cds.start==1], reftable$cds.id[reftable$cds.end==reftable$length]) # List of complete CDS
    ol_start = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="start", select="all"))[,1] # Genes overlapping contig starts
    if (any(ol_start)) {
        if (!is.na(reftable[ol_start,"5utr.start"])){
            reftable[ol_start,"5utr.start"] = reftable[ol_start,"5utr.start"] + 3 # Truncate the first 3 bases
        } else {
            reftable[ol_start,"chr.coding.start"] = reftable[ol_start,"chr.coding.start"] + 3 # Truncate the first 3 bases
            reftable[ol_start,"cds.start"] = reftable[ol_start,"cds.start"] + 3 # Truncate the first 3 bases
        }
    }

    ol_end = as.data.frame(GenomicRanges::findOverlaps(transc_gr, chrs_gr, type="end", select="all"))[,1] # Genes overlapping contig starts
    if (any(ol_end)) {
        if (!is.na(reftable[ol_start,"3utr.end"])){
            reftable[ol_start,"3utr.end"] = reftable[ol_start,"3utr.end"] - 3 # Truncate the first 3 bases
        } else {
            reftable[ol_end,"chr.coding.end"] = reftable[ol_end,"chr.coding.end"] - 3 # Truncate the first 3 bases
            reftable[ol_end,"cds.end"] = reftable[ol_end,"cds.end"] - 3 # Truncate the first 3 bases
        }
    }

    if (any(c(ol_start,ol_end))) {
        warning(sprintf("The following genes were found to start or end at the first or last base of their contig. Since dndscv needs trinucleotide contexts for all coding bases, codons overlapping ends of contigs have been trimmed. Affected genes: %s.", paste(reftable[unique(c(ol_start,ol_end)),"gene.name"], collapse=", ")))
    }

    # Selecting the longest complete CDS for every gene (required when there are multiple alternative transcripts per unique gene name)
    
    cds_table = unique(reftable[,c(1:3,9)])
    cds_table = cds_table[order(cds_table$gene.name, -cds_table$length), ] # Sorting CDS from longest to shortest
    cds_table = cds_table[(cds_table$length %% 3)==0, ] # Removing CDS of length not multiple of 3
    cds_table = cds_table[cds_table$cds.id %in% fullcds, ] # Complete CDS
    reftable = reftable[reftable$cds.id %in% fullcds, ] # Complete CDS
    gene_list = unique(cds_table$gene.name)
    
    reftable = reftable[order(reftable$chr, reftable$chr.coding.start), ]
    cds_split = split(reftable, f=reftable$cds.id)
    gene_split = split(cds_table, f=cds_table$gene.name)

    ## 2. Building the RefCDS object
    message("[2/3] Building the RefCDS object...")
    
    # Subfunction to extract the coding sequence
    get_CDSseq = function(gr, strand) {
        cdsseq = strsplit(paste(as.vector(Rsamtools::scanFa(genomefile, gr)),collapse=""),"")[[1]]
        if (strand==-1) {
            cdsseq = rev(seqinr::comp(cdsseq,forceToLower=F))
        }
        return(cdsseq)
    }
    
    # Subfunction to extract essential splice site sequences
    # Definition of essential splice sites: (5' splice site: +1,+2,+5; 3' splice site: -1,-2)
    get_splicesites = function(cds) {
        splpos = numeric(0)
        if (nrow(cds)>1) { # If the CDS has more than one exon
            if (cds[1,10]==1) { # + strand
                spl5prime = cds[-nrow(cds),6] # Exon end before splice site
                spl3prime = cds[-1,5] # Exon start after splice site
                splpos = unique(sort(c(spl5prime+1, spl5prime+2, spl5prime+5, spl3prime-1, spl3prime-2)))
            } else if (cds[1,10]==-1) { # - strand
                spl5prime = cds[-1,5] # Exon end before splice site
                spl3prime = cds[-nrow(cds),6] # Exon start after splice site
                splpos = unique(sort(c(spl5prime-1, spl5prime-2, spl5prime-5, spl3prime+1, spl3prime+2)))
            }
        }
        return(splpos)
    }
    
    # Subfunction to extract the essential splice site sequence
    get_spliceseq = function(gr, strand) {
        spliceseq = unname(as.vector(Rsamtools::scanFa(genomefile, gr)))
        if (strand==-1) {
            spliceseq = seqinr::comp(spliceseq,forceToLower=F)
        }
        return(spliceseq)
    }

    # Subfunction to define 5'UTR intervals
    get_5UTRposition = function(strand) {
        if (strand==1 & !is.na(cds[1, 12])) { # + strand and if 5'UTR start position exists
            utr5pos = seq(from=cds[1, 12], to=cds[1,13])
        }
        if (strand==-1 & !is.na(cds[nrow(cds),13])) {
            utr5pos = seq(from=cds[nrow(cds), 12], to=cds[nrow(cds),13])
        } else {
            utr5pos = numeric(0)
        }
        return(utr5pos)
    }

    # Subfunction to define 5'UTR intervals
    get_3UTRposition = function(strand) {
        if (strand==1 & !is.na(cds[nrow(cds), 14])) { # + strand and if 5'UTR start position exists
            utr3pos = seq(from=cds[nrow(cds), 14], to=cds[nrow(cds), 15])
        }
        if (strand==-1 & !is.na(cds[1, 15])) {
            utr3pos = seq(from=cds[1, 14], to=cds[1, 15])
        } else {
            utr3pos = numeric(0)
        }
        return(utr3pos)
    }

    # Subfunction to extract 5'UTR sequences
    get_5UTRseq = function(strand, shift=0) { # shift: number of nts to shift back - (to 5') or forward + (to 3') (e.g. +3 shift 3 nts forward)
        if (strand==1 & !is.na(cds[1, 12])) { # + strand and if 5'UTR start position exists
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[1, 12]+shift, cds[1,13]+shift)) # range of 5'UTR nts
            UTRseq = unname(as.vector(Rsamtools::scanFa(genomefile, gr))) 
        }
        if (strand==-1 & !is.na(cds[nrow(cds),13])) { # - strand and if 5'UTR start position exists
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[nrow(cds), 12]+shift, cds[nrow(cds),13]+shift)) # range of 5'UTR nts
            UTRseq = strsplit(paste(as.vector(Rsamtools::scanFa(genomefile, gr)),collapse=""),"")[[1]]
            UTRseq = seqinr::comp(UTRseq,forceToLower=F)
        }
        return(UTRseq)
    }

    # Subfunction to extract 3'UTR sequences
    get_3UTRseq = function(strand, shift=0) { # shift: number of nts to shift back - (to 5') or forward + (to 3') (e.g. +3 shift 3 nts forward)
        if (strand==1 & !is.na(cds[nrow(cds), 14])) { # + strand and if 3'UTR start position exists
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[nrow(cds), 14]+shift, cds[nrow(cds)+shift,15])) # range of 3'UTR nts
            UTRseq = unname(as.vector(Rsamtools::scanFa(genomefile, gr)))
        }
        if (strand==-1 & !is.na(cds[1, 15])) { # - strand and if 3'UTR start position exists
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[1, 14]+shift, cds[1,15]+shift))
            UTRseq = strsplit(paste(as.vector(Rsamtools::scanFa(genomefile, gr)),collapse=""),"")[[1]] # range of 3'UTR nts
            UTRseq = seqinr::comp(UTRseq,forceToLower=F)
        }
        return(UTRseq)
    }

    # Initialising and populating the RefCDS object
    
    RefCDS = array(list(NULL), length(gene_split)) # Initialising empty object
    invalid_genes = rep(0, length(gene_split)) # Initialising empty object
    
    for (j in 1:length(gene_split)) {
        
        gene_cdss = gene_split[[j]]
        h = keeptrying = 1
        
        while (h<=nrow(gene_cdss) & keeptrying) {
            
            pid = gene_cdss[h,3]
            cds = cds_split[[pid]]
            strand = cds[1,10]
            chr = cds[1,4]
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5], cds[,6]))
            cdsseq = get_CDSseq(gr,strand)
            pseq = seqinr::translate(cdsseq, numcode = numcode)
            if (all(pseq[-length(pseq)]!="*") & all(cdsseq!="N")) { # A valid CDS has been found (no stop codons inside the CDS excluding the last codon) and no "N" nucleotides
                
                # Essential splice sites
                splpos = get_splicesites(cds) # Essential splice sites
                if (length(splpos)>0) { # CDSs with a single exon do not have splice sites
                    gr_spl = GenomicRanges::GRanges(chr, IRanges::IRanges(splpos, splpos))
                    splseq = get_spliceseq(gr_spl, strand)
                }

                # 5'UTR and 3'UTR 
                utr5pos = get_5UTRposition(strand); utr3pos = get_3UTRposition(strand)
                if (length(utr5pos)>0) { # CDSs with 5'UTR
                    utr5seq = get_5UTRseq(strand, shift=0)
                }
                if (length(utr3pos)>0) { # CDSs with 3'UTR
                    utr3seq = get_3UTRseq(strand, shift=0)
                }
                
                # Obtaining the splicing sequences and the coding and splicing sequence contexts
                if (strand==1) {
                    
                    cdsseq1up = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]-1, cds[,6]-1)), strand)
                    cdsseq1down = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]+1, cds[,6]+1)), strand)

                    if (length(splpos)>0) {
                        splseq1up = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos-1, splpos-1)), strand)
                        splseq1down = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos+1, splpos+1)), strand)
                    }

                    if (length(utr5pos)>0) { 
                        utr5seq1up = get_5UTRseq(strand, shift=-1)
                        utr5seq1down = get_5UTRseq(strand, shift=1)
                    }

                    if (length(utr3pos)>0) { 
                        utr3seq1up = get_3UTRseq(strand, shift=-1)
                        utr3seq1down = get_3UTRseq(strand, shift=1)
                    }
                    
                } else if (strand==-1) {
                    
                    cdsseq1up = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]+1, cds[,6]+1)), strand)
                    cdsseq1down = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]-1, cds[,6]-1)), strand)

                    if (length(splpos)>0) {
                        splseq1up = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos+1, splpos+1)), strand)
                        splseq1down = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos-1, splpos-1)), strand)
                    }

                    if (length(utr5pos)>0) { 
                        utr5seq1up = get_5UTRseq(strand, shift=1)
                        utr5seq1down = get_5UTRseq(strand, shift=-1)
                    }

                    if (length(utr3pos)>0) { 
                        utr3seq1up = get_3UTRseq(strand, shift=1)
                        utr3seq1down = get_3UTRseq(strand, shift=-1)
                    }
                    
                }
                
                # Annotating the CDS in the RefCDS database
                
                RefCDS[[j]]$gene_name = gene_cdss[h,2]
                RefCDS[[j]]$gene_id = gene_cdss[h,1]
                RefCDS[[j]]$protein_id = gene_cdss[h,3]
                RefCDS[[j]]$CDS_length = gene_cdss[h,4]
                RefCDS[[j]]$chr = cds[1,4]
                RefCDS[[j]]$strand = strand
                RefCDS[[j]]$intervals_cds = unname(as.matrix(cds[,5:6]))
                RefCDS[[j]]$intervals_splice = splpos
                RefCDS[[j]]$intervals_5utr = utr5pos
                RefCDS[[j]]$intervals_3utr = utr3pos
                
                RefCDS[[j]]$seq_cds = Biostrings::DNAString(paste(cdsseq, collapse=""))
                RefCDS[[j]]$seq_cds1up = Biostrings::DNAString(paste(cdsseq1up, collapse=""))
                RefCDS[[j]]$seq_cds1down = Biostrings::DNAString(paste(cdsseq1down, collapse=""))

                if (length(utr5pos) > 0) {
                    RefCDS[[j]]$utr5seq = Biostrings::DNAString(paste(utr5seq, collapse=""))
                    RefCDS[[j]]$utr5seq1up = Biostrings::DNAString(paste(utr5seq1up, collapse=""))
                    RefCDS[[j]]$utr5seq1down = Biostrings::DNAString(paste(utr5seq1down, collapse=""))
                }
                
                if (length(utr3pos) > 0) {
                    RefCDS[[j]]$utr3seq = Biostrings::DNAString(paste(utr3seq, collapse=""))
                    RefCDS[[j]]$utr3seq1up = Biostrings::DNAString(paste(utr3seq1up, collapse=""))
                    RefCDS[[j]]$utr3seq1down = Biostrings::DNAString(paste(utr3seq1down, collapse=""))
                }
                
                if (length(splpos)>0) { # If there are splice sites in the gene
                    RefCDS[[j]]$seq_splice = Biostrings::DNAString(paste(splseq, collapse=""))
                    RefCDS[[j]]$seq_splice1up = Biostrings::DNAString(paste(splseq1up, collapse=""))
                    RefCDS[[j]]$seq_splice1down = Biostrings::DNAString(paste(splseq1down, collapse=""))
                }
                
                keeptrying = 0 # Stopping the while loop
            }
            h = h+1
        }
        if (keeptrying) {
            invalid_genes[j] = 1 # No valid CDS was found for this gene and the gene will be removed from the RefCDS object
        }
        if (round(j/1000)==(j/1000)) { message(sprintf('    %0.3g%% ...', round(j/length(gene_split),2)*100)) }
    }
    
    RefCDS = RefCDS[!invalid_genes] # Removing genes without a valid CDS

    ## 3. L matrices: number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context
    message("[3/3] Calculating the impact of all possible coding changes...")
    
    nt = c("A","C","G","T")
    trinuc_list = paste(rep(nt,each=16,times=1), rep(nt,each=4,times=4), rep(nt,each=1,times=16), sep="")
    trinuc_ind = structure(1:64, names=trinuc_list)

    trinuc_subs = NULL; for (j in 1:length(trinuc_list)) { trinuc_subs = c(trinuc_subs, paste(trinuc_list[j], paste(substr(trinuc_list[j],1,1), setdiff(nt,substr(trinuc_list[j],2,2)), substr(trinuc_list[j],3,3), sep=""), sep=">")) }
    trinuc_subsind = structure(1:192, names=trinuc_subs)

    # Precalculating a 64x64 matrix with the functional impact of each codon transition (1=Synonymous, 2=Missense, 3=Nonsense)
    impact_matrix = array(NA, dim=c(64,64))
    colnames(impact_matrix) = rownames(impact_matrix) = trinuc_list
    for (j in 1:64) {
        for (h in 1:64) {
            from_aa = seqinr::translate(strsplit(trinuc_list[j],"")[[1]], numcode = numcode)
            to_aa = seqinr::translate(strsplit(trinuc_list[h],"")[[1]], numcode = numcode)
            # Annotating the impact of the mutation
            if (to_aa == from_aa){ 
                impact_matrix[j,h] = 1
            } else if (to_aa == "*"){
                impact_matrix[j,h] = 3
            } else if ((to_aa != "*") & (from_aa != "*") & (to_aa != from_aa)){
                impact_matrix[j,h] = 2
            } else if (from_aa=="*") {
                impact_matrix[j,h] = NA
            }
        }
    }
    
    for (j in 1:length(RefCDS)) {
        
        L = array(0, dim=c(192,6))
        cdsseq = as.character(as.vector(RefCDS[[j]]$seq_cds))
        cdsseq1up = as.character(as.vector(RefCDS[[j]]$seq_cds1up))
        cdsseq1down = as.character(as.vector(RefCDS[[j]]$seq_cds1down))
        
        # 1. Exonic mutations
        
        ind = rep(1:length(cdsseq), each=3)
        old_trinuc = paste(cdsseq1up[ind], cdsseq[ind], cdsseq1down[ind], sep="")
        new_base = c(sapply(cdsseq, function(x) nt[nt!=x]))
        new_trinuc = paste(cdsseq1up[ind], new_base, cdsseq1down[ind], sep="")
        codon_start = rep(seq(1,length(cdsseq),by=3),each=9)
        old_codon = paste(cdsseq[codon_start], cdsseq[codon_start+1], cdsseq[codon_start+2], sep="")
        pos_in_codon = rep(rep(1:3, each=3), length.out=length(old_codon))
        aux = strsplit(old_codon,"")
        new_codon = sapply(1:length(old_codon), function(x) { new_codonx = aux[[x]]; new_codonx[pos_in_codon[x]] = new_base[x]; return(new_codonx) } )
        new_codon = paste(new_codon[1,], new_codon[2,], new_codon[3,], sep="")
        
        imp = impact_matrix[(trinuc_ind[new_codon]-1)*64 + trinuc_ind[old_codon]]
        matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
        
        # Synonymous
        matrix_ind = table(matrind[which(imp==1)])
        L[as.numeric(names(matrix_ind)), 1] = matrix_ind
        
        # Missense
        matrix_ind = table(matrind[which(imp==2)])
        L[as.numeric(names(matrix_ind)), 2] = matrix_ind
        
        # Nonsense
        matrix_ind = table(matrind[which(imp==3)])
        L[as.numeric(names(matrix_ind)), 3] = matrix_ind
        
        # 2. Splice site mutations
        if (length(RefCDS[[j]]$intervals_splice)>0) {
            splseq = as.character(as.vector(RefCDS[[j]]$seq_splice))
            splseq1up = as.character(as.vector(RefCDS[[j]]$seq_splice1up))
            splseq1down = as.character(as.vector(RefCDS[[j]]$seq_splice1down))
            old_trinuc = rep(paste(splseq1up, splseq, splseq1down, sep=""), each=3)
            new_trinuc = paste(rep(splseq1up, each=3), c(sapply(splseq, function(x) nt[nt!=x])), rep(splseq1down,each=3), sep="")
            matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
            matrix_ind = table(matrind)
            L[as.numeric(names(matrix_ind)), 4] = matrix_ind
        }

        if (length(RefCDS[[j]]$intervals_5utr)>0) {
            utr5seq = as.character(as.vector(RefCDS[[j]]$utr5seq))
            utr5seq1up = as.character(as.vector(RefCDS[[j]]$utr5seq1up))
            utr5seq1down = as.character(as.vector(RefCDS[[j]]$utr5seq1down))
            old_trinuc = rep(paste(utr5seq1up, utr5seq, utr5seq1down, sep=""), each=3)
            new_trinuc = paste(rep(utr5seq1up, each=3), c(sapply(utr5seq, function(x) nt[nt!=x])), rep(utr5seq1down,each=3), sep="")
            matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
            matrix_ind = table(matrind)
            L[as.numeric(names(matrix_ind)), 5] = matrix_ind
        }

        if (length(RefCDS[[j]]$intervals_3utr)>0) {
            utr3seq = as.character(as.vector(RefCDS[[j]]$utr3seq))
            utr3seq1up = as.character(as.vector(RefCDS[[j]]$utr3seq1up))
            utr3seq1down = as.character(as.vector(RefCDS[[j]]$utr3seq1down))
            old_trinuc = rep(paste(utr3seq1up, utr3seq, utr3seq1down, sep=""), each=3)
            new_trinuc = paste(rep(utr3seq1up, each=3), c(sapply(utr3seq, function(x) nt[nt!=x])), rep(utr3seq1down,each=3), sep="")
            matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
            matrix_ind = table(matrind)
            L[as.numeric(names(matrix_ind)), 6] = matrix_ind
        }

        RefCDS[[j]]$L = L # Saving the L matrix
        if (round(j/1000)==(j/1000)) { message(sprintf('    %0.3g%% ...', round(j/length(gene_split),2)*100)) }
    }

    ## Saving the reference GenomicRanges object

    aux = unlist(sapply(1:length(RefCDS), function(x) t(cbind(x,rbind(RefCDS[[x]]$intervals_cds,cbind(RefCDS[[x]]$intervals_splice,RefCDS[[x]]$intervals_splice),
    cbind(RefCDS[[x]]$intervals_5utr,RefCDS[[x]]$intervals_5utr),cbind(RefCDS[[x]]$intervals_3utr,RefCDS[[x]]$intervals_3utr))))))
    df_genes = as.data.frame(t(array(aux,dim=c(3,length(aux)/3))))
    colnames(df_genes) = c("ind","start","end")
    df_genes$chr = unlist(sapply(1:length(RefCDS), function(x) rep(RefCDS[[x]]$chr,nrow(RefCDS[[x]]$intervals_cds)+length(RefCDS[[x]]$intervals_splice)+length(RefCDS[[x]]$intervals_5utr)+length(RefCDS[[x]]$intervals_3utr))))
    df_genes$gene = sapply(RefCDS, function(x) x$gene_name)[df_genes$ind]
    
    gr_genes = GenomicRanges::GRanges(df_genes$chr, IRanges::IRanges(df_genes$start, df_genes$end))
    GenomicRanges::mcols(gr_genes)$names = df_genes$gene
    
    save(RefCDS, gr_genes, file=outfile)

}

# EOF



    

    

