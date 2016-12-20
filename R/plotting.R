#' Create box plot of normalized read counts
#'
#' @param result result object of panelcn.mops
#' @param sampleName name of the test sample that should be displayed
#' @param countWindows data.frame with contents of a BED file as returned by
#' getWindows
#' @param selectedGenes vector of names of genes of interest that should be
#' displayed or NULL if all genes are of interest. Default = NULL
#' @param showGene integer indicating which of the genes of interest to plot
#' @param showLegend flag to indicate whether to display a legend with the
#' names of the test samples. Default = TRUE
#' @param exonRange vector of 2 positive integers to limit box plot to a
#' certain range of exons or NULL
#' @return generates a boxplot of the normalized read counts
#' @examples
#' data(panelcn.mops)
#' sampleNames <- colnames(elementMetadata(test))
#' selectedGenes <- "ATM"
#' plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1], 
#'            countWindows = countWindows, selectedGenes = selectedGenes, 
#'            showGene = 1)
#' @export
plotBoxplot <- function(result, sampleName, countWindows, selectedGenes = NULL,
                        showGene = 1, showLegend = TRUE, exonRange = NULL) {

    if (missing(countWindows)) {
        stop("\"countWindows\" need to be specified.")
    }
    
    if (missing(result)) {
        stop("\"result\" needs to be specified.")
    }
    
    if (missing(sampleName)) {
        stop("\"sampleName\" needs to be specified.")
    }
    
    if (is.null(selectedGenes)) {
        message("All genes selected.")
        selectedGenes <- unique(countWindows$gene)
    }
    
    if (!is.null(showGene) && length(showGene == 1) &&
        showGene <= length(selectedGenes)) {
        geneidx <- showGene
    } else {
        geneidx <- 1
        message(paste0("Setting showGene=", showGene,
                        " not possible - using gene ", selectedGenes[geneidx]))
    }
    if (sampleName != result@sampleNames[1]) {
        message(paste0("sampleName: ", sampleName, " not equal to name ",
                        "specified in result which is: ", 
                        result@sampleNames[1]))
    }
    
    selectedGenes <- selectedGenes[geneidx]
    geneWindows <- countWindows[which(countWindows$gene %in% selectedGenes),]
    geneWindowsPaste <- paste(geneWindows$chromosome, geneWindows$start,
                                geneWindows$end, sep="_")
    plotData <- as.matrix(result@normalizedData)[geneWindowsPaste,]
    geneWindowsData <- geneWindows[which(geneWindowsPaste %in%
                                            rownames(plotData)),]
    genes <- geneWindowsData$gene
    exons <- sapply(strsplit(geneWindowsData[,4], "[.]"), "[[", 2)


    startLabels <- paste(exons, " (", 1:length(exons), ")")

    if (is.null(exonRange) || length(exonRange) != 2 || 0 %in% exonRange) {
        exonRange <- c(1, length(startLabels))
    }

    m <- length(sampleName)
    n <- ncol(plotData)

    ylim <- c(0, max(plotData[exonRange[1]:exonRange[2],]))

    if (showLegend) {
        #par(mar=c(5.1, 4.1, 4.1, 16.1), xpd=NA)
        #par(xpd=T, mar=par()$mar+c(0,0,length(filelist),0))
        ylim <- ylim*1.05
    }

    boxplot(t(plotData[exonRange[1]:exonRange[2],]), ylim=ylim, xaxt="n",
            ylab="normalized read counts", bty='L', main=unlist(selectedGenes))

    axis(1, at=1:(exonRange[2]-exonRange[1]+1),
            labels=startLabels[exonRange[1]:exonRange[2]], las=2)
    col_vec <- c(rainbow(m), rep("black", n-m))
    for (j in n:1) {
        points(plotData[exonRange[1]:exonRange[2],j], col=col_vec[j], pch=19,
                cex=0.5)
    }

    if (showLegend) {
        #"topright", inset=c(-1,0)
        #0, max(plotData[exonRange[1]:exonRange[2],])+30+length(filelist)
        legend("topright", sampleName, pch = 19, col=col_vec[1:m], cex=1)
    }

    #par(mar=c(5, 4, 4, 2) + 0.1)
}
