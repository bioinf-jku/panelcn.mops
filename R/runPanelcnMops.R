#' @title Full copy number detection for targeted NGS panel data for 
#' multiple samples
#' @description This function performs first quality control and runs 
#' panelcn.mops for CNV detection on all test samples.
#' @param XandCB GRanges object of combined  read counts of test samples and
#' control samples as returned by countBamListInGRanges
#' @param testiv vector of indices of test samples in XandCB. Default = c(1)
#' @param countWindows data.frame with contents of a BED file as returned by
#' getWindows
#' @param selectedGenes vector of names of genes of interest or NULL if all 
#' genes are of interest. Default = NULL
#' @param I vector of positive real values containing the expected fold change
#' of the copy number classes. Length of this vector must be equal to the
#' length of the "classes" parameter vector. For targeted NGS panel data
#' the default is c(0.025,0.57,1,1.46,2)
#' @param normType type of the normalization technique. Each samples'
#' read counts are scaled such that the total number of reads are comparable
#' across samples. Options are "mean","median","poisson", "quant", and "mode"
#' Default = "quant"
#' @param sizeFactor parameter for calculating the size factors for
#' normalization. Options are "mean","median", "quant", and "mode".
#' Default = "quant"
#' @param qu Quantile of the normType if normType is set to "quant".
#' Real value between 0 and 1. Default = 0.25
#' @param quSizeFactor Quantile of the sizeFactor if sizeFactor is set to
#' "quant". 0.75 corresponds to "upper quartile normalization".
#' Real value between 0 and 1. Default = 0.75
#' @param norm the normalization strategy to be used. If set to 0 the read
#' counts are not normalized and cn.mops does not model different coverages.
#' If set to 1 the read counts are normalized. If set to 2 the read counts are
#' not normalized and cn.mops models different coverages. Default = 1.
#' @param priorImpact positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will be
#' assumed to have copy number 2. Default = 1
#' @param minMedianRC segments with median read counts over
#' all samples < minMedianRC are excluded from the analysis
#' @param maxControls integer reflecting the maximal numbers of controls to 
#' use. If set to 0 all highly correlated controls are used. Default = 25
#' @param corrThresh threshold for selecting highly correlated controls. 
#' Default = 0.99
#' @param sex either "mixed", "male", or "female" reflecting the sex of
#' all samples (test and control)
#' @return list of instances of "CNVDetectionResult"
#' @import S4Vectors
#' @importClassesFrom cn.mops CNVDetectionResult
#' @examples
#' data(panelcn.mops)
#' XandCB <- test
#' elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
#'                                     elementMetadata(control))
#' resultlist <- runPanelcnMops(XandCB, countWindows = countWindows)
#' @export
runPanelcnMops <- function(XandCB, testiv = c(1), countWindows, 
                            selectedGenes = NULL, 
                            I = c(0.025, 0.57, 1, 1.46, 2), 
                            normType = "quant", sizeFactor = "quant", 
                            qu = 0.25, quSizeFactor = 0.75,
                            norm = 1, priorImpact = 1, minMedianRC = 30, 
                            maxControls = 25, corrThresh = 0.99,
                            sex = "mixed") {

    if (missing(countWindows)) {
        stop("\"countWindows\" need to be specified.")
    }
    if(!(sex %in% c("mixed", "male", "female"))) {
        message(paste0("Setting sex=", sex, " not possible - ",
                        "using sex=\"mixed\""))
    }

    if (is.null(selectedGenes)) {
        message("All genes selected.")
        selectedGenes <- c()
    }

    XandCB@elementMetadata <- XandCB@elementMetadata[,c(testiv, 
                (1:ncol(XandCB@elementMetadata))[-testiv])]
    testiv <- 1:length(testiv)
    
    sampleNames <- colnames(XandCB@elementMetadata)
    message(paste0("Analyzing sample(s) ", sampleNames[testiv], "\n"))

    XandCBMatrix <- as.matrix(XandCB@elementMetadata)

    ## quality control
    maxRC <- apply(XandCBMatrix, 1, max)
    medianRC <- apply(XandCBMatrix, 1, median)
    sampleMedian <- apply(XandCBMatrix, 2, median)
    sampleThresh <- median(sampleMedian[-testiv])*0.55
#    sampleThresh <- mean(sampleMedian[-testiv]) - 2*sd(sampleMedian[-testiv])
    message(paste("new sampleThresh", sampleThresh))
    poorQual <- which(medianRC < minMedianRC)
    highRC <- which(maxRC >= 5000 & maxRC < 25000)
    veryHighRC <- which(maxRC >= 25000)
    poorSamples <- which(sampleMedian < sampleThresh)

    for (h in highRC) {
        for (s in seq_len(ncol(XandCBMatrix))) {
            XandCB@elementMetadata[h,s] <- XandCBMatrix[h,s]/10
        }
    }

    for (h in veryHighRC) {
        for (s in seq_len(ncol(XandCBMatrix))) {
            XandCB@elementMetadata[h,s] <- XandCBMatrix[h,s]/100
        }
    }

    colnames(XandCB@elementMetadata) <- sampleNames
    
    if (length(highRC) > 0){
        message(paste0("Had to reduce read counts for exon ",
                        countWindows[highRC,]$name,"\n"))
    }
    if (length(veryHighRC) > 0){
        message(paste0("Had to reduce read counts for exon ",
                        countWindows[veryHighRC,]$name,"\n"))
    }


    if (length(poorQual) > 0) {
        message(paste("Cannot use exon", countWindows[poorQual,]$name, "\n"))
    }

    XChr <- c(which(countWindows$chromosome=="chrX" |
                    countWindows$chromosome=="X"))

    if (length(XChr) > 0) {
        if (sex=="mixed") {
            message(paste0("Ignoring X-chromosomal exons ",
                            "(sex is mixed/unknown).\n"))
        } else {
            message(paste0("All females or all males selected. ", 
                            "Chromosome X treated like autosomes."))
            XChr <- c()
        }
        if (sex=="male") {
            message("Male: Note that CN2 is actually CN1 for chromosome X.")
        }
    }

    YChr <- c(which(countWindows$chromosome=="chrY" |
                    countWindows$chromosome=="Y"))

    if (length(YChr) > 0) {
        message(paste0("Ignoring Y-chromosomal exons."))
    }

    ignoreExons <- unique(c(poorQual, XChr, YChr))
    subsetIdx <- rep(TRUE, nrow(countWindows))
    subsetIdx[ignoreExons] <- FALSE
    usedExons <- seq_len(nrow(countWindows))[-ignoreExons]
    if (length(ignoreExons) > 0) {
        countWindows <- countWindows[-ignoreExons,]
    }
    countWindows <- countWindows[order(suppressWarnings(
                        as.numeric(countWindows[,1])), countWindows[,2]),]
    
    if (length(selectedGenes) > 0) {
        geneInd <- c()
        for (g in selectedGenes) {
            geneIndTemp <- which(countWindows$gene==g)
            if (length(geneIndTemp) == 0) {
                message(paste0("Gene ", g, " not in \"countWindows\""))
            }
            geneInd <- c(geneInd, geneIndTemp)
        }
        
        
        if (length(geneInd) == 0) {
            stop(paste0("At least one of the \"selectedGenes\" needs to be ", 
                        "in \"countWindows\"."))
        }
    } else {
        geneInd <- NULL
    }
    

    poorDBSamples <- poorSamples[!(poorSamples %in% testiv)]
    poorTestSamples <- poorSamples[poorSamples %in% testiv]


    if (length(poorSamples) > 0) {
        message(paste("Ignoring bad control sample", sampleNames[poorDBSamples],
                        "\n"))
    }

    if (length(poorTestSamples) > 0) {
        message(paste("Bad test sample", sampleNames[poorTestSamples], "\n"))
    }
    poorSamples <- poorDBSamples

    if (length(poorSamples) > 0) {
        XandCB <- XandCB[,-poorSamples]
        sampleNames <- sampleNames[-poorSamples]
        colnames(XandCB@elementMetadata) <- sampleNames
    }

    ii <- 1
    resultlist <- list()
    for (t in testiv) {
        message(paste0("\nAnalyzing sample ", sampleNames[t], "\n"))
        controli <- seq_len(ncol(XandCB@elementMetadata))[-testiv]
        dup <- grep(sampleNames[t], sampleNames[-testiv])

        if (length(dup) > 0) {
            message("Removing test sample from control samples\n")
            controli <- controli[-dup]
        }
        result <- panelcn.mops(subset(XandCB[,c(t,controli)], subsetIdx),
                            testi = 1, geneInd = geneInd, I = I, 
                            priorImpact = priorImpact,
                            normType = normType, sizeFactor = sizeFactor, 
                            qu = qu, quSizeFactor = quSizeFactor, norm = norm,
                            maxControls = maxControls, corrThresh = corrThresh)
        
        resultlist[[ii]] <- result
        ii <- ii + 1
    }

    return(resultlist)
}


#' Test data included in panelcn.mops
#' @name test
#' @docType data
#' @title GRanges object of countWindows with read counts for a test sample as 
#' elementMetadata. 
#' @description The object was created using the function 
#' countBamListInGRanges with the enclosed countWindows object, a subset of a 
#' BAM file provided by the 1000 Genomes Project and the read.width parameter 
#' set to 150.
#' @keywords data
#' @examples
#' data(panelcn.mops)
#' test
#' @author Gundula Povysil
NULL

#' Control data included in panelcn.mops
#' @name control
#' @docType data
#' @title GRanges object of countWindows with read counts for control samples 
#' as elementMetadata. 
#' @description The object was created using the function 
#' countBamListInGRanges with the enclosed countWindows object, a subset of 
#' BAM files provided by the 1000 Genomes Project and the read.width parameter 
#' set to 150.
#' @keywords data
#' @examples
#' data(panelcn.mops)
#' control
#' @author Gundula Povysil
NULL

#' Data included in panelcn.mops
#' @name countWindows
#' @docType data
#' @title result object of getWindows - a data.frame with the contents of 
#' the provided BED file with an additional gene name and exon name column
#' @examples
#' data(panelcn.mops)
#' countWindows
#' @keywords data
#' @author Gundula Povysil
NULL

#' Result data included in panelcn.mops
#' @name resultlist
#' @docType data
#' @title result object of runPanelcnMops - a list of instances of 
#' "CNVDetectionResult"
#' @keywords data
#' @examples
#' data(panelcn.mops)
#' resultlist
#' @author Gundula Povysil
NULL

#' Data included in panelcn.mops
#' @name read.width
#' @docType data
#' @title read width used for calculating RCs of test and control
#' @keywords data
#' @examples
#' data(panelcn.mops)
#' read.width
#' @author Gundula Povysil
NULL

