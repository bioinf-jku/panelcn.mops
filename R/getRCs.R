#' Convert BED file into data.frame of count windows
#'
#' @param filename filename of the BED file with absolute or relative path
#' (structure of BED file without header: chromosome, exon start, exon end,
#' exon name)
#' @param chr indicates whether naming contains chr prefix
#' @return a data.frame with the contents of the BED file with an additional
#' gene name and exon name column
#' @examples
#' bed <- list.files(system.file("extdata", package = "panelcn.mops"),
#'                         pattern = ".bed$", full.names = TRUE)
#' countWindows <- getWindows(bed)
#' @importFrom utils read.csv
#' @export
getWindows <- function(filename, chr = FALSE) {
    data <- read.csv(filename, sep="\t", header = FALSE)
    if (ncol(data) < 4 || is.na(data[1,2])) {
        stop("BED file needs to have gene name in 4th column and no header.")
    }
    if (!is.numeric(data$V2) || !is.numeric(data$V3)) {
        stop(paste0("2nd and 3rd column of BED file need to be numeric. ",
                "Make sure all columns are separated by tabs."))
    }
    data <- data[,1:4]
    data$V4 <- as.character(data$V4)
    data$V1 <- as.character(data$V1)
    data <- data[!duplicated(paste(data$V1, data$V2, data$V3, sep="_")),]
    firstColumnAsList <- strsplit(data[,4],"[.]")
    firstColumnAsVector <- sapply(firstColumnAsList,'[',1)
    data[,5] <- firstColumnAsVector
    secondColumnAsVector <- sapply(firstColumnAsList,'[',2)
    if (length(grep("E", secondColumnAsVector))) {
        data[,6] <- secondColumnAsVector
    } else {
        data[,6] <- NA
    }

    names(data) <- c("chromosome", "start", "end", "name", "gene", "exon")
    if (!chr && any(grepl("chr", data$chromosome))) {
        data[,1] <- sub("chr", "", data[,1])
        message(paste0("naming without chr prefix chosen, but BED contains ",
                        "chr -> removing chr"))
    } else if (chr && !any(grepl("chr", data$chromosome))) {
        data[,1] <- paste0("chr", data[,1])
        message(paste0("naming with chr prefix chosen, but BED does not",
                        " contain chr -> adding chr"))
    }
    return(data)
}


#' Get read counts for a list of BAM files and given count windows
#'
#' @param bam.files list with absolute or relative paths to BAM files
#' @param countWindows data.frame with contents of a BED file as returned by
#' getWindows
#' @param read.width read.width parameter for countBamInGRanges or FALSE if
#' actual read width should be extracted from BAM file
#' @param ... additional parameters
#' @return a GRanges object over the countWindows with read counts for each
#' sample as elementMetadata
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBam
#' @examples
#' bed <- list.files(system.file("extdata", package = "panelcn.mops"),
#'                     pattern = ".bed$", full.names = TRUE)
#' countWindows <- getWindows(bed)
#' testbam <- list.files(system.file("extdata", package = "panelcn.mops"),
#' pattern = ".bam$", full.names = TRUE)
#' test <- countBamListInGRanges(countWindows = countWindows,
#'                                 bam.files = testbam, read.width = 150)

countBamListInGRanges <- function(bam.files, countWindows, read.width = 150, 
                                    ...){
    RCs <- list()
    GR <- GRanges(countWindows[,1], IRanges(countWindows[,2], 
                                            countWindows[,3]))

    for (i in 1:length(bam.files)) {
        message(paste("Processing", basename(bam.files[i]), "...", i, "/",
                        length(bam.files)))
        indexed <- (file.exists(paste(bam.files[i],".bai",sep="")) | 
                    file.exists(gsub(".bam",".bai", bam.files[i])))
        
        if (indexed){
            if (read.width == FALSE) {
                message("get.width")
                RCs[[i]] <- .countBamInGRanges(bam.file=bam.files[i], GR ,
                                                min.mapq=1, get.width = TRUE,
                                                ...)
            } else {
                RCs[[i]] <- .countBamInGRanges(bam.file=bam.files[i], GR ,
                                                min.mapq=1, 
                                                read.width = read.width, ...)
            }
        } else {
            message(paste0("No index file (.bai) found for BAM file", 
                            bam.files[i], ". RCs cannot be calculated"))
            return()
        }
    }

    message("finished processing samples")
    M <- matrix(unlist(RCs), nrow=nrow(countWindows))

    if (all(M == 0)) {
        message(paste0("All read counts are 0. Please make sure ",
                        "that you have the right BAM file and BED file."))
    }

#    colnames(M) <- basename(bam.files)
#    print(colnames(M))
    RC <- GR
    RC@elementMetadata <- DataFrame(M)
    colnames(RC@elementMetadata) <- basename(bam.files)
    return(RC)
}

#' countBamInGRanges from package exomeCopy
.countBamInGRanges <- function(bam.file, granges, min.mapq = 1, read.width = 1,
                                stranded.start = FALSE,
                                get.width = FALSE, remove.dup = FALSE) {
    rds.counts <- integer(length(granges))
    seq.names <- unique(as.character(GenomicRanges::seqnames(granges)))
    seq.names.in.bam <- names(Rsamtools::scanBamHeader(bam.file)[[1]]$targets)

    if ((sum(grepl(pattern = '^chr', seq.names)) > 0) &&
        (sum(grepl(pattern = '^chr', seq.names.in.bam)) == 0)) {
        warning('countWindows contain chr prefix to the chromosome name but
                BAM does not -> removing chr.')
        seq.names <- substr(seq.names, 4, nchar(seq.names))
        GenomeInfoDb::seqlevels(granges) <- 
            substr(GenomeInfoDb::seqlevels(granges), 4, 
                    nchar(GenomeInfoDb::seqlevels(granges)))
    } else if ((sum(grepl(pattern = '^chr', seq.names.in.bam)) > 0) &&
        (sum(grepl(pattern = '^chr', seq.names)) == 0)) {
        warning('BAM contains chr prefix to the chromosome name but
                countWindows do not -> adding chr.')
        seq.names <- paste0("chr", seq.names)
        GenomeInfoDb::seqlevels(granges) <- 
            paste0("chr", GenomeInfoDb::seqlevels(granges))
    }

    for (seq.name in seq.names) {
        if (seq.name %in% seq.names.in.bam) {
            granges.subset <-
                granges[GenomicRanges::seqnames(granges) == seq.name]
            strand(granges.subset) <- "*"
            scan.what <- c("pos", "mapq")
            if (stranded.start | get.width) {
                scan.what <- c(scan.what, "qwidth")
            }
            if (stranded.start | remove.dup) {
                scan.what <- c(scan.what, "strand")
            }
            rds <-
                Rsamtools::scanBam(bam.file,
                    param = Rsamtools::ScanBamParam(what = scan.what,
                                                which = range(granges.subset)))
            if (min.mapq > 0) {
                mapq.test <- rds[[1]]$mapq >= min.mapq & !is.na(rds[[1]]$mapq)
            }
            else {
                mapq.test <- rep(TRUE, length(rds[[1]]$mapq))
            }
            if (remove.dup) {
                if (get.width) {
                    mapq.test <- mapq.test &
                        !duplicated(rds[[1]]$pos +
                                    as.numeric(rds[[1]]$strand)/10 +
                                                rds[[1]]$qwidth/1e+05)
                }
                else {
                    mapq.test <- mapq.test & !duplicated(rds[[1]]$pos +
                                                as.numeric(rds[[1]]$strand)/10)
                }
            }
            if (sum(mapq.test) > 0) {
                if (stranded.start) {
                    rds.ranges <-
                        GenomicRanges::GRanges(seq.name,
                            IRanges::IRanges(start =
                                ifelse(rds[[1]]$strand[mapq.test] %in%
                                    c("+", "*"), rds[[1]]$pos[mapq.test],
                                    rds[[1]]$pos[mapq.test] +
                                    rds[[1]]$qwidth[mapq.test] - read.width),
                                end = ifelse(rds[[1]]$strand[mapq.test] %in%
                                    c("+", "*"), rds[[1]]$pos[mapq.test] +
                                    read.width - 1, rds[[1]]$pos[mapq.test] +
                                    rds[[1]]$qwidth[mapq.test] - 1)))
                }
                else if (get.width) {
                    rds.ranges <-
                        GenomicRanges::GRanges(seq.name,
                            IRanges::IRanges(start = rds[[1]]$pos[mapq.test],
                                width = rds[[1]]$qwidth[mapq.test]))
                }
                else {
                    rds.ranges <- GenomicRanges::GRanges(seq.name,
                            IRanges::IRanges(start = rds[[1]]$pos[mapq.test],
                                                width = read.width))
                }
                rds.counts.seq.name <-
                    GenomicRanges::countOverlaps(granges.subset, rds.ranges)
                rds.counts[as.logical(GenomicRanges::seqnames(granges) ==
                                            seq.name)] <- rds.counts.seq.name
            }
            else {
                rds.counts[as.logical(GenomicRanges::seqnames(granges) ==
                                            seq.name)] <- 0
            }
        }
        else {
            rds.counts[as.logical(GenomicRanges::seqnames(granges) ==
                                        seq.name)] <- 0
        }
    }
    if (sum(rds.counts) == 0) {
        message("No reads found with minimum mapping quality")
    }
    rds.counts
}


#' Split (larger) ROIs into multiple smaller (overlapping) bins and create 
#' new BED file
#' 
#' @param oldBedFile filename of the BED file with absolute or relative path
#' (structure of BED file without header: chromosome, exon start, exon end,
#' exon name)
#' @param newBedFile filename of the new BED file that should be created
#' @param limit ROIs larger than limit will be split
#' @param bin size of bins (in bp) the ROIs will be split into
#' @param shift no. of bp between start positions of adjacent bins 
#' @param chr indicates whether naming contains chr prefix
#' @return generates a new BED file with (larger) ROIs split into smaller bins
#' @importFrom utils write.table
#' @export
#' @examples
#' bed <- list.files(system.file("extdata", package = "panelcn.mops"),
#'                     pattern = ".bed$", full.names = TRUE)
#' splitROIs(bed, "newBed.bed")


splitROIs <- function(oldBedFile, newBedFile, limit = 0, bin = 100, shift = 50, 
                        chr = FALSE) {
  
    if (!(is.numeric(limit) & limit >= 0 & length(limit)==1)) {
        stop("\"limit\" must be numeric, larger or equal 0 and of length 1.")
    }
    
    if (!(is.numeric(bin) & bin > 0 & length(bin)==1)) {
        stop("\"bin\" must be numeric, larger than 0 and of length 1.")
    }
    
    if (!(is.numeric(shift) & shift > 0 & length(shift)==1)) {
        stop("\"shift\" must be numeric, larger than 0 and of length 1.")
    }
    
    message(paste("reading from bedFile", basename(oldBedFile)))

    windows <- getWindows(oldBedFile, chr)
    
    
    message(paste("No. of original ROIs:", nrow(windows)))
    

    starts <- c()
    ends <- c()
    chrom <- c()
    exon <- c()
    geneName <- c()

    for (i in 1:nrow(windows)) {
    cw_row <- windows[i,]
    if ((cw_row$end - cw_row$start) > limit & 
        (cw_row$end - cw_row$start) > bin) {
        starts <- c(starts, seq(cw_row$start, (cw_row$end - bin + 1), shift))
        ends <- c(ends, seq((cw_row$start + bin - 1), cw_row$end, shift))
        
        diff <- ends[length(ends)] - cw_row$end
        if (diff != 0) {
            stopifnot(diff < 0)
            ends[length(ends)] <- cw_row$end
        }
        nrows <- length(ends) - length(chrom)
        chrom <- c(chrom, rep(cw_row$chrom, nrows))
        geneName <- c(geneName, rep(cw_row$gene, nrows))
        exon <- c(exon, rep(cw_row$exon, nrows))
    } else {
        starts <- c(starts, cw_row$start)
        ends <- c(ends, cw_row$end)
        
        chrom <- c(chrom, cw_row$chrom)
        geneName <- c(geneName, cw_row$gene)
        exon <- c(exon, cw_row$exon)
    }
    }

    if (sum(grepl(pattern = '^chr', chrom)) == 0) {
        chromn <- paste0("chr", chrom)
    } else {
        chromn <- chrom
    }

    geneName <- paste(geneName, exon, chromn, starts, ends, sep = ".")
    message(paste("Created", length(chrom), "ROIs"))

    windows <- data.frame(chrom = chrom, start = starts, end = ends, 
                            name = geneName, stringsAsFactors = FALSE)

    write.table(windows, file = newBedFile, sep = "\t", quote = FALSE, 
                col.names = FALSE, row.names = FALSE)

}


