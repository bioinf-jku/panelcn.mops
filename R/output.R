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
        tempTable$sampleName <- as.character(tempTable$sampleName)
        tempTable <- tempTable[which(tempTable$sampleName == 
                                        tempTable$sampleName[1]),]
        tmpNames <- c(paste("X", sampleNames, sep = ""), sampleNames)
        tmpRows <- which(tempTable$sampleName %in% tmpNames)
        
        if (length(tmpRows) > 0) {
            tempTable <- tempTable[tmpRows,]
            message(paste(unique(tempTable$sampleName), "\n"))
            
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
            exonNr <- rep(0, nrow(tempTable))
            
            ## low Qual
            lowQ <- rep("", nrow(tempTable))
            badexi <- result@params$badROI
            badexn <- row.names(result@normalizedData)[badexi]
            
            ##build initial table
            message("Building table...")
            tempTable <- data.frame(tempTable, genes, exonNr, exons, RC, 
                                    medianRC, RC.norm, medianRC.norm, lowQ, 
                                    stringsAsFactors = FALSE)
            
            ## nrExons
            nrExons <- nrow(tempTable) / length(unique(tempTable$sampleName))
            
            ## select by gene
            geneWindows <- countWindows[which(countWindows$gene %in% 
                                                selectedGenes),]
            
            tempTable <- tempTable[which(paste(tempTable$seqnames,
                                                tempTable$start, 
                                                tempTable$end, sep = "_") %in%
                                        paste(geneWindows$chromosome, 
                                            geneWindows$start, geneWindows$end,
                                                sep = "_")),]
            
            ## used in function
            medianRC <- apply(as.matrix(XandCB@elementMetadata), 1, median)
            medianRCNorm <- apply(as.matrix(result@normalizedData), 1, median)
            
            currSample <- tempTable[1,]$sampleName
            message(currSample)
            
            tempWindows <- paste(tempTable$seqnames,
                                    tempTable$start, 
                                    tempTable$end, sep = "_")
            
            gw.rows <- geneWindows[which(paste(geneWindows$chromosome,
                                                geneWindows$start, 
                                                geneWindows$end, sep = "_") %in%
                                                tempWindows),]
            row.names(gw.rows) <- paste(gw.rows$chromosome,
                                        gw.rows$start, 
                                        gw.rows$end, sep = "_")
            # reorder
            gw.rows <- gw.rows[tempWindows,]
            
            tempTable$genes <- gw.rows$gene
            tempTable$exons <- gw.rows$name
            
            
            XandCBred <- XandCB[which(paste(as.vector(seqnames(XandCB)), 
                                            start(XandCB), end(XandCB), 
                                            sep = "_") %in% tempWindows),
                                which(colnames(XandCB@elementMetadata) 
                                        == currSample)]
            
            RCs <- XandCBred@elementMetadata[,1]
            redNames <- paste(as.vector(seqnames(XandCBred)), start(XandCBred),
                                end(XandCBred), sep = "_")
            names(RCs) <- redNames
            # make sure that RCs have right order
            tempTable$RC <- RCs[tempWindows]
            
            tempMedianRC <- medianRC[which(paste(as.vector(seqnames(XandCB)), 
                                                start(XandCB), end(XandCB), 
                                                sep = "_") %in% tempWindows)]
            names(tempMedianRC) <- redNames
            
            tempTable$medianRC <- tempMedianRC[tempWindows]
            
            tempTable$RC.norm <- 
                round(result@normalizedData[which(
                    row.names(result@normalizedData) 
                                %in% tempWindows),
                    which(colnames(result@normalizedData) == 
                            currSample)])[tempWindows]
            
            
            tempTable$medianRC.norm <- 
                round(medianRCNorm[which(row.names(result@normalizedData) 
                                            %in% tempWindows)])[tempWindows]
            
            
            tempTable$exonNr <- gw.rows$exon
            
            
            ccn <- as.character(cn[tempWindows, currSample])
            tempTable$CN <- ccn
            
            badw <- which(tempWindows %in% badexn)
            if (length(badw) > 0) {
                tempTable$lowQ[badw] <- "lowQual"
            }
            
            
            # GRanges adds X to sampleNames that start with number
            Xidx <- which(!(tempTable$sampleName %in% sampleNames))
            tempTable[Xidx,]$sampleName <-
                substr(tempTable[Xidx,]$sampleName, start = 2,
                        stop = nchar(as.character(tempTable[Xidx,]$sampleName)))
            
            tempTable <- data.frame("Sample" = tempTable$sampleName,
                                "Chr" = tempTable$seqnames, 
                                "Gene" = tempTable$genes,
                                "Exon" = tempTable$exons,
                                # "Exon" = tempTable$exonNr,
                                "Start" = tempTable$start, 
                                "End" = tempTable$end,
                                "RC" = tempTable$RC, 
                                "medRC" = tempTable$medianRC,
                                "RC.norm" = tempTable$RC.norm,
                                "medRC.norm" = tempTable$medianRC.norm,
                                "lowQual" = tempTable$lowQ,
                                "CN" = tempTable$CN)
            
            
            resulttable[[samp]] <- tempTable
            
        } else {
            message(paste0("Sample ", tempTable$sampleName[1], 
                            " not selected."))
        }
        message("Finished")



    }
    return(resulttable)
}
