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
#' @param ylimup numeric, maximum RC is multiplied by this value to calculate 
#' second value of ylim. Default = 1.15
#' @param thresh numeric threshold for plotting fold change areas 
#' E.g. thresh = 0.4 plots a green rectangle above (1 + 0.4)*median for each 
#' boxplot and a red rectangle below (1 - 0.4)*median. Default of zero does 
#' not plot any colored areas.
#' @return generates a boxplot of the normalized read counts
#' @examples
#' data(panelcn.mops)
#' sampleNames <- colnames(elementMetadata(test))
#' selectedGenes <- "ATM"
#' plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1], 
#'             countWindows = countWindows, selectedGenes = selectedGenes, 
#'             showGene = 1)
#' @importFrom graphics points axis boxplot legend par rect
#' @importFrom grDevices rainbow rgb
#' @export
plotBoxplot <- function(result, sampleName, countWindows, selectedGenes = NULL,
                        showGene = 1, showLegend = TRUE, exonRange = NULL, 
                        ylimup = 1.15, thresh = 0) {

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
    
    if (!is.numeric(thresh)) {
        stop("\"thresh\" needs to be numeric")
    }
    
    if (!is.numeric(ylimup)) {
        stop("\"ylimup\" needs to be numeric")
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
    countWindowsPaste <- paste(countWindows$chromosome, countWindows$start,
                               countWindows$end, sep="_")
    plotData <- as.matrix(result@normalizedData)
    countWindowsNew <- countWindows[which(countWindowsPaste %in%
                                             rownames(plotData)),]
    
    # dummy ROI necessary for genes with only 1 ROI
    dummy <- seq_len(nrow(countWindowsNew))[-which(countWindowsNew$gene %in% selectedGenes)][1]
    if (is.na(dummy)) {
        stop("Gene not in result - nothing to plot!")
    }
    geneWindows <- countWindowsNew[c(dummy, which(countWindowsNew$gene %in% selectedGenes)),]
    geneWindowsPaste <- paste(geneWindows$chromosome, geneWindows$start,
                                geneWindows$end, sep="_")
    plotData <- plotData[geneWindowsPaste,]
    
    genes <- geneWindows$gene[-1]
    exons <- sapply(strsplit(geneWindows[,4], "[.]"), "[[", 2)[-1]


    startLabels <- paste(exons, " (", seq_along(exons), ")")

    if (is.null(exonRange) || length(exonRange) != 2 || 0 %in% exonRange) {
        exonRange <- c(1, length(startLabels))
    }

    par(mar=c(6, 4, 4, 2) + 0.1)

    m <- length(sampleName)
    n <- ncol(plotData)

    ylim <- c(0, (max(plotData[(exonRange[1]+1):(exonRange[2]+1),]))*ylimup)

    if (abs(exonRange[2]-exonRange[1]) > 1) {
        boxplot(t(plotData[(exonRange[1]+1):(exonRange[2]+1),]), ylim=ylim, 
                xaxt="n", ylab="normalized read counts", bty='L', 
                main=unlist(selectedGenes))
    } else {
        boxplot(plotData[(exonRange[1]+1):(exonRange[2]+1),], ylim=ylim, 
                xaxt="n", ylab="normalized read counts", bty='L', 
                main=unlist(selectedGenes))
    }
    
    if (thresh) {
        x0s <- seq_along((exonRange[1]):(exonRange[2])) - 0.4
        x1s <- seq_along((exonRange[1]):(exonRange[2])) + 0.4
        y0s <- apply(plotData[(exonRange[1]+1):(exonRange[2]+1),], 1, 
                        median)*(1 + thresh)
        y1s <- apply(plotData[(exonRange[1]+1):(exonRange[2]+1),], 1, 
                        median)*(1 - thresh)

#        segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "green", lty=2)
#        segments(x0 = x0s, x1 = x1s, y0 = y1s, col = "red", lty=2)
        rect(xleft = x0s, xright = x1s, ybottom = y0s, ytop = ylim[2], 
                col = rgb(0, 255, 0, 50, maxColorValue=255), lty = "blank")
        rect(xleft = x0s, xright = x1s, ybottom = ylim[1], ytop = y1s, 
                col = rgb(255, 0, 0, 50, maxColorValue=255), lty = "blank")
    }

    axis(1, at=seq_len(abs(exonRange[2]-exonRange[1])+1),
            labels=startLabels[(exonRange[1]):(exonRange[2])], las=2)
    col_vec <- c(rainbow(m), rep("black", n-m))
    for (j in n:1) {
        points(plotData[(exonRange[1]+1):(exonRange[2]+1),j], col=col_vec[j], 
                pch=19, cex=0.5)
    }

    if (showLegend) {
        #"topright", inset=c(-1,0)
        #0, max(plotData[exonRange[1]:exonRange[2],])+30+length(filelist)
        legend("topright", sampleName, pch = 19, col=col_vec[seq_len(m)], cex=1)
    }

    par(mar=c(5, 4, 4, 2) + 0.1)
}
