#' Creates a user readable result table for the test samples of the
#' genes of interest
#'
#' @param resultlist result object of runPanelcnMops
#' @param XandCB GRanges object of combined  read counts of test samples and
#' control samples as returned by getRCRanges or countBamListInGRanges
#' @param countWindows data.frame with contents of a BED file as returned by
#' getWindows
#' @param selectedGenes vector of names of genes of interest that should be
#' displayed or NULL if all genes are of interest. Default = NULL
#' @param sampleNames names of the test samples (basename of the BAM files)
#' @return a data.frame containing the results for the test samples within the
#' genes of interest
#' @examples
#' data(panelcn.mops)
#' XandCB <- test
#' elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
#'                                 elementMetadata(control))
#' sampleNames <- colnames(elementMetadata(test))
#' selectedGenes <- "ATM"
#' resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, 
#'                                     countWindows = countWindows, 
#'                                     selectedGenes = selectedGenes, 
#'                                     sampleNames = sampleNames)
#' @export
createResultTable <- function(resultlist, XandCB, countWindows,
                                selectedGenes = NULL, sampleNames){

    if (missing(countWindows)) {
        stop("\"countWindows\" need to be specified.")
    }
    
    if (missing(resultlist)) {
        stop("\"resultlist\" needs to be specified.")
    }
    
    if (missing(XandCB)) {
        stop("\"XandCB\" needs to be specified.")
    }
    
    if (missing(sampleNames)) {
        stop("\"sampleNames\" need to be specified.")
    }
    
    message(paste0("Calculating results for sample(s) ", sampleNames, '\n'))
    
    if (is.null(selectedGenes)) {
        message("All genes selected.")
        selectedGenes <- unique(countWindows$gene)
    }

    resulttable <- list()
    
    for (samp in 1:length(resultlist)) {
        result <- resultlist[[samp]]
    
        tempTable <- as.data.frame(cn.mops::cnvs(result))

        ## CNs
        cn <- result@integerCopyNumber
        
        ## seqnames
        tempTable$seqnames <- as.character(tempTable$seqnames)
       
        ## sampleNames
        # message("before sample name selection...")
        tempTable$sampleName <- as.character(tempTable$sampleName)
        # GRanges adds X to sampleNames that start with number
        tmpNames <- c(paste("X",sampleNames, sep=""), sampleNames)
        tmpRows <- which(tempTable$sampleName %in% tmpNames)

        tempTable <- tempTable[tmpRows,]
        
        tmpCols <- which(colnames(cn) %in% tmpNames)
        names <- colnames(cn)[tmpCols]
        cn <- as.data.frame(cn[,tmpCols])
        colnames(cn) <- names
        
        ## RC and median RC
        # message("before rc...")
        RC <- rep(NA, nrow(tempTable))
        medianRC <- rep(NA, nrow(tempTable))
        RC.norm <- rep(NA, nrow(tempTable))
        medianRC.norm <- rep(NA, nrow(tempTable))

        ## genes and exons
        genes <- rep(NA, nrow(tempTable))
        exons <- rep(NA, nrow(tempTable))
        #exonsPerGenes <- aggregate(name ~ gene, data=countWindows, length)
        exonNr <- rep(0, nrow(tempTable))

        ## low Qual
        lowQ <- rep(NA, nrow(tempTable))
        badexi <- result@params$badROI
        badexn <- row.names(result@normalizedData)[badexi]
        
        ##build initial table
        message("Building table...")
        tempTable <- data.frame(tempTable, genes, exonNr, exons, RC, medianRC,
                                RC.norm, medianRC.norm, lowQ)

        ##get starts in order of normalized matrix -- order is different
        normDataStarts <- as.numeric(sapply(strsplit(
                            row.names(result@normalizedData), "[_]"), "[[", 2))
        normDataEnds <- as.numeric(sapply(strsplit(
            row.names(result@normalizedData), "[_]"), "[[", 3))
        

        ## nrExons
        nrExons <- nrow(tempTable) / length(unique(tempTable$sampleName))

        ## select by gene
        geneWindows <- countWindows[which(countWindows$gene %in% 
                                                selectedGenes),]
        tempTable <- tempTable[which(paste(tempTable$seqnames,
                                            tempTable$start, 
                                            tempTable$end, sep="_") %in%
                                    paste(geneWindows$chromosome, 
                                            geneWindows$start, geneWindows$end,
                                            sep="_")),]

        ## used in function
        medianRC <- apply(as.matrix(XandCB@elementMetadata), 1, median)
        medianRCNorm <- apply(as.matrix(result@normalizedData), 1, median)

        for (i in 1:nrow(tempTable)) {
            chr <- tempTable[i,]$seqnames
            start <- tempTable[i,]$start
            end <- tempTable[i,]$end
            startEnd <- paste(start, end, sep="_")
            
            currSample <- tempTable[i,]$sampleName

            gw.row <- geneWindows[which(paste(geneWindows$chromosome,
                                            geneWindows$start, 
                                            geneWindows$end, sep="_") ==
                                        paste(chr, start, end, sep="_")),]

            tempTable[i,]$genes <- gw.row[,5]
            tempTable[i,]$exons <- gw.row[,4]

            tempTable[i,]$RC <- XandCB[which(paste(start(XandCB), end(XandCB), 
                                                    sep="_") == startEnd),
                                        which(colnames(XandCB@elementMetadata) 
                                            == currSample)]@elementMetadata[1,1]
            tempTable[i,]$medianRC <- medianRC[which(paste(start(XandCB), 
                                                        end(XandCB), sep="_") 
                                                    == startEnd)]

            tempTable[i,]$RC.norm <-
                round(result@normalizedData[which(paste(normDataStarts, 
                                                        normDataEnds, sep="_") 
                                                    == startEnd),
                            which(colnames(result@normalizedData)==currSample)])
            tempTable[i,]$medianRC.norm <-
                round(medianRCNorm[which(paste(normDataStarts, normDataEnds, 
                                                sep="_") == startEnd)])

            tempTable[i,]$exonNr <- gw.row$exon
            
            postSeg <- paste(gw.row$chromosome, gw.row$start, gw.row$end, 
                                sep="_")
            ccn <- as.character(cn[postSeg, currSample])
            tempTable[i,]$CN <- ccn
            
            if (postSeg %in% badexn) {
                tempTable[i,]$lowQ <- "lowQual"
            } else {
                tempTable[i,]$lowQ <- ""
            }
        }
        # to have order from temptable, not from countWindows
        tempGenes <- unique(tempTable$genes)
        tempSamples <- unique(tempTable$sampleName)

        # GRanges adds X to sampleNames that start with number
        Xidx <- which(!(tempTable$sampleName %in% sampleNames))
        tempTable[Xidx,]$sampleName <-
            substr(tempTable[Xidx,]$sampleName, start = 2,
                    stop = nchar(as.character(tempTable[Xidx,]$sampleName)))

        tempTable <- data.frame("Sample"=tempTable$sampleName,
                                "Chr"=tempTable$seqnames, 
                                "Gene"=tempTable$genes,
                                "Exon"=tempTable$exons,
                                # "Exon"=tempTable$exonNr,
                                "Start"=tempTable$start, "End"=tempTable$end,
                                "RC"=tempTable$RC, "medRC"=tempTable$medianRC,
                                "RC.norm"=tempTable$RC.norm,
                                "medRC.norm"=tempTable$medianRC.norm,
                                "lowQual"=tempTable$lowQ,
                                "CN"=tempTable$CN)

        message("Finished")

        resulttable[[samp]] <- tempTable
    }
    return(resulttable)
}
