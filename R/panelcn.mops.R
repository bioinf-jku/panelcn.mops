#' main function of the package panelcn.mops
#'
#' @param input either an instance of "GRanges" or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An entry is
#' the read count of a sample in the genomic region.
#' @param testi positive integer that gives the index of the test sample in 
#' input. Default = 1
#' @param geneInd vector of indices of rows input that are within target genes. 
#' These regions are not considered for chosing correlated reference samples.
#' If NULL, all regions are considered for the correlation. Default = NULL
#' @param classes vector of characters of the same length as the parameter
#' vector "I". One vector element must be named "CN2". The names reflect the
#' labels of the copy number classes.
#' Default = c("CN0","CN1","CN2","CN3","CN4").
#' @param I vector of positive real values containing the expected fold change
#' of the copy number classes. Length of this vector must be equal to the
#' length of the "classes" parameter vector. For human copy number polymorphisms
#' the default is c(0.025,0.5,1,1.5,2).
#' @param priorImpact positive real value that reflects how strong the prior
#' assumption affects the result. The higher the value the more samples will be
#' assumed to have copy number 2. Default = 1.
#' @param cyc positive integer that sets the number of cycles for the algorithm.
#' Usually after less than 15 cycles convergence is reached. Default = 20.
#' @param normType type of the normalization technique. Each samples'
#' read counts are scaled such that the total number of reads are comparable
#' across samples. Options are "mean", "median", "poisson", "quant", and "mode".
#' Default = "quant".
#' @param sizeFactor parameter for calculating the size factors for
#' normalization. Options are "mean", "median", "quant", and "mode".
#' Default = "quant".
#' @param qu Quantile of the normType if normType is set to "quant".
#' Real value between 0 and 1. Default = 0.25.
#' @param quSizeFactor Quantile of the sizeFactor if sizeFactor is set to
#' "quant". 0.75 corresponds to "upper quartile normalization".
#' Real value between 0 and 1. Default = 0.75.
#' @param norm the normalization strategy to be used. If set to 0 the read
#' counts are not normalized and cn.mops does not model different coverages.
#' If set to 1 the read counts are normalized. If set to 2 the read counts are
#' not normalized and cn.mops models different coverages. Default = 1.
#' @param minReadCount if all samples are below this value the algorithm will
#' return the prior knowledge. This prevents that the algorithm from being
#' applied to segments with very low coverage. Default = 5.
#' @param maxControls integer reflecting the maximal numbers of controls to 
#' use. If set to 0 all highly correlated controls are used. Default = 25
#' @param useMedian flag indicating whether "median" instead of "mean" of a
#' segment should be used for the CNV call. Default = FALSE.
#' @param returnPosterior flag that decides whether the posterior probabilities
#' should be returned. The posterior probabilities have a dimension of samples
#' times copy number states times genomic regions and therefore consume a lot of
#' memory. Default = TRUE.
#' @return an instance of "CNVDetectionResult".
#' @importFrom cn.mops normalizeChromosomes
#' @importFrom GenomeInfoDb sortSeqlevels
#' @importFrom methods new
#' @importFrom graphics boxplot
#' @examples
#' data(panelcn.mops)
#' XandCB <- test
#' elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
#'                                     elementMetadata(control))
#' result <- panelcn.mops(XandCB)

#' @export

panelcn.mops <- function(input, testi = 1, geneInd=NULL,
                            classes = c("CN0", "CN1", "CN2", "CN3", "CN4"),
                            I = c(0.025, 0.5, 1, 1.5, 2), priorImpact = 1, 
                            cyc = 20,
                            normType = "quant", sizeFactor = "quant", qu = 0.25,
                            quSizeFactor = 0.75, norm = 1, minReadCount = 5,
                            maxControls = 25,
                            useMedian = FALSE, returnPosterior = TRUE){
    # CHECK INPUT

    if(class(input)=="GRanges"){
        inputType <- "GRanges"
        input <- sortSeqlevels(input)
        input <- GenomicRanges::sort(input)
        X <- IRanges::as.matrix(IRanges::values(input))
        if (ncol(X)==1){
            stop("It is not possible to run cn.mops on only ONE sample.\n")
        }
        if (length(IRanges::unique(strand(input))) >1){
            stop(paste("Different strands found in GRanges object. Please make",
                                "read counts independent of strand."))
        }
        chr <- as.character(seqnames(input))
        start <- start(input)
        end <- end(input)
        rownames(X) <- paste(chr,start,end,sep="_")

        irAllRegions <- IRanges(start,end)
        grAllRegions <- GRanges(chr,irAllRegions,seqinfo=seqinfo(input))
        grAllRegions <- GenomeInfoDb::sortSeqlevels(grAllRegions)
        names(irAllRegions) <- NULL
    } else if (is.matrix(input)){
        if (nrow(input)> 1){
            inputType <- "DataMatrix"
            X <- input
            X <- matrix(as.numeric(X),nrow=nrow(X))
            colnames(X) <- colnames(input)
            chr <- rep("undef",nrow(X))
            irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
            grAllRegions <- GRanges(chr,irAllRegions)
        } else{
            inputType <- "DataMatrix"
            chr <- "undef"
            irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
            grAllRegions <- GRanges(chr,irAllRegions)
            parallel <- 0
        }
    } else if (is.vector(input)) {
        inputType <- "DataMatrix"
        X <- matrix(input,nrow=1)
        X <- matrix(as.numeric(X),nrow=nrow(X))
        chr <- "undef"
        irAllRegions <- IRanges(start=1:nrow(X),end=1:nrow(X))
        grAllRegions <- GRanges(chr,irAllRegions)
        parallel <- 0
    }else{
        stop("GRanges object or read count matrix needed as input.")
    }
    

    if (!(is.numeric(testi)) & testi < 0 & length(testi)!=1)
        stop("\"testi\" must be numeric, larger than 0 and of length 1.")
    else if (testi > ncol(X))
        stop("\"testi\" must be smaller than the total number of samples.")
    else if (testi != 1) {
        message("Rearranging X because test sample is not in column 1")
        testi <- as.integer(testi)
        X <- X[,c(testi, 1:(testi-1), testi+(1:ncol(X)))]
    }

    if (!is.numeric(geneInd) & !is.null(geneInd)) 
        stop("\"geneInd\" must be numeric or NULL.")
    
    if (any(X<0) | any(!is.finite(X))){
        stop("All values must be greater or equal zero and finite.\n")
    }
    if (!is.numeric(I)) stop("\"I\" must be numeric.")
    if (!is.character(classes)) stop("\"classes\" must be character.")
    if (length(I)!=length(classes)){
        stop("I and classes must have same length!")
    }
    if (!("CN2" %in% classes)){stop("One element of classes must be CN2 .\n")}
    if (!(is.numeric(priorImpact) & length(priorImpact)==1))
        stop("\"priorImpact\" be must numeric and of length 1.")
    if (!(is.numeric(cyc) & length(cyc)==1))
        stop("\"cyc\" must be numeric and of length 1.")
    #if (!(is.numeric(parallel) & length(parallel)==1))
    #    stop("\"parallel\" must be numeric and of length 1.")
    if (!(normType %in% c("mean","min","median","quant","mode","poisson"))){
        stop(paste("Set normalization to \"mean\"",
                    "\"min\", \"median\", \"quant\" or \"mode\"."))
    }
    if (!(is.numeric(qu) & length(qu)==1))
        stop("\"qu\" must be numeric and of length 1.")
    if (!(is.numeric(quSizeFactor) & length(quSizeFactor)==1))
        stop("\"quSizeFactor\" must be numeric and of length 1.")
    if (is.logical(norm))
        norm <- as.integer(norm)
    if (!(norm %in% c(0,1,2)))
        stop("\"norm\" must be 0,1 or 2.")
    #if (!(is.numeric(minWidth) & length(minWidth)==1))
    #    stop("\"minWidth\" must be numeric and of length 1.")
    #if (!is.character(segAlgorithm)){
    #    stop("\"segAlgorithm\" must be \"fastseg\" or \"DNAcopy\"!")
    #}
    
    if (!(is.numeric(maxControls) & maxControls > 0 & length(maxControls)==1)) {
        stop("\"maxControls\" must be numeric, larger than 0 and of length 1.")
    } else {
        message(paste0("\"maxControls\" is set to ", maxControls))
    }
    
    if (!(is.numeric(minReadCount) & length(minReadCount)==1))
        stop("\"minReadCount\" must be numeric and of length 1.")
    if (!is.logical(returnPosterior))
        stop("\"returnPosterior\" must be logical.")

    ########################################################################

    params <- list("panelcn.mops", I, classes,
                    priorImpact, cyc,
                    sizeFactor, normType, qu, quSizeFactor,
                    minReadCount, useMedian, "CN2", version,
                    paste("..."))
    names(params) <- c("method", "I", "classes",
                        "priorImpact", "cyc",
                        "sizeFactor", "normType", "qu", "quSizeFactor",
                        "minReadCount", "useMedian", "mainClass",
                        "version", "SegmentationParams")

    ########################################################################
    nROIs <- nrow(X)
    nSamples <- ncol(X)
    n <- length(I)
    
    sampleNames <- colnames(X)
    if(is.null(sampleNames)) {
        sampleNames <- as.character(1:nSamples)
    }
    sampleNamesC <- sampleNames[-testi]
    sampleNamesT <- sampleNames[testi]
    
    # WARNINGS ####################################################

    if (nROIs < 100){
        warning(paste("Normalization might not be applicable",
                            "for this small number of segments."))
    }

    # selecting controls ###################################

    goodROI <- 1:nROIs
    goodROI <- setdiff(goodROI, geneInd)
    if (length(goodROI)>0) {
        Xgood <- X[goodROI,]
    } else {
        Xgood <- X
    }
    
    K <- cor(Xgood)
    newO <- order(K[,1], decreasing = TRUE)
    K <- K[newO,]
    
    message(paste0(K[,1], "\n"))
    
    conti <- which(K[,1] > 0.99)

    conti <- setdiff(conti, 1)
    if (length(conti) < 8) {
        message("Low correlation.\n")
        if (nSamples > 8) {
            conti <- order(K[-c(1),1], decreasing = TRUE)[1:8]+1
        } else {
            conti <- 2:nSamples
        }
    }
    if (maxControls > 0) {
        if (length(conti) > maxControls) {
            conti <- conti[1:maxControls]
        }
    }
    message(paste0("Selected ", length(conti), " out of ", 
                    length(sampleNamesC), " controls:\n"))
    message(paste0(sampleNames[newO[conti]], "\n"))
    
    sampleNames <- c(sampleNamesT, sampleNames[newO[conti]])
    X <- X[,sampleNames]
    nSamples <- ncol(X)
    
    
    # NORMALIZING ################################################


    if (norm==0){
        X.norm <- X
        cov <- rep(1,nSamples)
    } else if (norm==1) {
        message("Normalizing...")
        X.norm <- normalizeChromosomes(X, normType=normType,
                                        sizeFactor=sizeFactor, qu=qu, 
                                        quSizeFactor=quSizeFactor)
        medianX.norm <- apply(as.matrix(X.norm), 1, median)
        ratios <- sapply(1:ncol(X.norm), function(i)
                            as.matrix(X.norm)[,i]/medianX.norm )
        b <- boxplot(ratios, plot=FALSE)
        c <- (b$stats[5,] - b$stats[1,])
        bad <- which(c > 0.5)

        badtest <- intersect(bad, testi)
        badcontrol <- setdiff(bad, testi)
        
        if (length(sampleNames) - length(badcontrol) < 5) {
            message("Too many bad control samples")
            cs <- sort(c)
            print(c)
            badcontrol <- c()
        }
        
        if (length(badcontrol) > 0) {
            message(paste("Removed bad control sample(s)",
                            sampleNames[badcontrol], "\n"))
            sampleNames <- sampleNames[-badcontrol]
            X <- X[,-badcontrol]
            message("Normalizing again...")
            X.norm <- normalizeChromosomes(X, normType=normType,
                                            sizeFactor=sizeFactor, qu=qu,
                                            quSizeFactor=quSizeFactor)
            message(paste("Remaining samples:",
                            sampleNames, "\n"))
            nSamples <- ncol(X)
        
            medianX.norm <- apply(as.matrix(X.norm), 1, median)
            ratios <- sapply(1:ncol(X.norm), function(i)
                as.matrix(X.norm)[,i]/medianX.norm )
            b <- boxplot(ratios, plot=FALSE)
            bad <- which((b$stats[5,] - b$stats[1,]) > 0.5)
            badtest <- intersect(bad, testi)
            badcontrol <- setdiff(bad, testi)
            if (length(badcontrol) > 0) {
                message(paste0("Still bad control"))
            }
        }
        if (length(badtest) > 0) {
            message(paste("Bad test sample", sampleNames[badtest], 
                            "\nTry using more controls!"))
        }
        
        be <- boxplot(t(X.norm), plot=FALSE)
        ce <- (be$stats[5,] - be$stats[1,])/be$stats[3,]
        badROI <- which(ce > 0.5)

        if (length(badROI) > 0) {
            message("Low quality exon(s):")
            message(paste0(rownames(X)[badROI], "\n"))
        }
        
        cov <- rep(1,nSamples)
    } else if (norm==2) {
        X.viz <- normalizeChromosomes(X, normType=normType,
                                        sizeFactor=sizeFactor, qu=qu, 
                                        quSizeFactor=quSizeFactor)
        X.norm <- X
        # robust estimates for the different coverages
        cov <- apply(X.norm,2,function(x) {
            mean(x[which(x>quantile(x,0.05) & x < quantile(x,0.95))])
        })
        if (median(cov)==0)
            stop("Median of the coverages is zero!")
        cov <- cov/median(cov)
    } else {
        stop("\"norm\" must be 0,1 or 2.")
    }
    params$cov <- cov

    
    chrOrder <- unique(chr) #unique chromosome names in alphabetic order
    chrBpts <- cumsum(table(chr)[chrOrder])
    # contains the idx where chromosomes start and end in X
    chrDf <- data.frame(start=c(1,chrBpts[-length(chrBpts)]+1),
                        end=chrBpts)
    rownames(chrDf) <- chrOrder
    
    
    # CN.MOPS COLUMNWISE

    # message("Starting local modeling, please be patient...    ")
    message("Computing results...")
    res <- list()

    res <- apply(X.norm[, , drop=FALSE], 1, cn.mops:::.cn.mopsC, I=I,
                    classes=classes, cov=cov, cyc=cyc, priorImpact=priorImpact,
                    minReadCount=minReadCount)

    ## Postprocess result
    L <- t(sapply(res,.subset2,1))
    rownames(L) <- rownames(X)
    colnames(L) <- classes
    params$L <- L
    A <- t(sapply(res,.subset2,2))
    rownames(A) <- rownames(X)
    colnames(A) <- classes
    params$A <- A
    CN <- t(sapply(res,.subset2,3))
    rownames(CN) <- rownames(X)
    colnames(CN) <- colnames(X)
    sINI <- t(sapply(res,.subset2,4))
    rownames(sINI) <- rownames(X)
    colnames(sINI) <- colnames(X)
    INI <- (sapply(res,.subset2,5))
    names(INI) <- rownames(X)
    
    if (returnPosterior){
        tt <- try(post <- array(dim=c(nROIs,n,nSamples)))
        if (inherits(tt,"try-error")){
            message("Posterior too large for array extent.")
            post <- array(NA,dim=c(1,1,1))
        } else {
            post.tmp <- t(lapply(res,.subset2,6))
            for (i in 1:nROIs){
                post[i, ,] <- post.tmp[[i]]
            }
            dimnames(post) <- list(rownames(X),classes,colnames(X))
            rm("post.tmp")
        }
    } else {
        post <- array(NA,dim=c(1,1,1))
    }
    rm("res")

    resSegmList <- list()
    segDf <- data.frame(stringsAsFactors=FALSE)

    callsS <- matrix(NA,nrow=nROIs,ncol=nSamples)
    colnames(callsS) <- colnames(X)

    for (chrom in chrOrder) {
        chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]


        resSegmList[[chrom]] <- apply(sINI[chrIdx, , drop=FALSE], 2,
                    function(x) {
                        nbr <- length(chrIdx)
                        data.frame("start"=1:nbr, "end"=1:nbr,
                                    "mean"=x[1:nbr], "median"=x[1:nbr])
                    })

        segDfTmp <- cbind(do.call(rbind,resSegmList[[chrom]]),
                            "sample"=rep(colnames(X),
                            sapply(resSegmList[[chrom]],nrow)))

        segDfTmp$chr <- chrom

        segDf <- rbind(segDf, segDfTmp)

        if (useMedian){
            callsS[chrIdx, ] <- matrix(rep(segDfTmp$median,
                                            segDfTmp$end - segDfTmp$start+1),
                                        ncol=nSamples)
        } else {
            callsS[chrIdx, ]<- matrix(rep(segDfTmp$mean,
                                            segDfTmp$end - segDfTmp$start+1),
                                        ncol=nSamples)
        }
    }

    segDf <- data.frame(segDf, "CN"=NA, stringsAsFactors=FALSE)

    #browser()

    colnames(segDf) <- c("start","end","mean","median","sample",
                                            "chr","CN")
    segDf <- segDf[ ,c("chr","start","end","sample","median","mean","CN")]


    segDfSubset <- segDf


    if (nrow(segDfSubset)>0){
        message("Creating CNVDetectionResult")

        # Assembly of result object
        r <- new("CNVDetectionResult")

        cnvrR <- reduce(GRanges(seqnames=segDfSubset$chr,
                                IRanges(segDfSubset$start,segDfSubset$end),
                                seqinfo=seqinfo(grAllRegions)))
        cnvrR <- GenomeInfoDb::sortSeqlevels(cnvrR)

        cnvrCN <- matrix(NA, ncol=nSamples, nrow=length(cnvrR))

        colnames(cnvrCN) <- colnames(X)

        sampleNames <- segDfSubset$sample

        if (inputType=="GRanges"){
            ir <- IRanges()
            irCNVR <- IRanges()
            for (chrom in chrOrder){
                #message(chrom)
                chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
                inputChr <- input[chrIdx]
                segDfSubsetChr <- subset(segDfSubset,chr==chrom)
                cnvrRChr <- cnvrR[which(as.character(
                    seqnames(cnvrR))==chrom)]
                if (nrow(segDfSubsetChr) >0){
                    ir <- c(ir,IRanges(start(inputChr)[
                        segDfSubsetChr$start],
                        end(inputChr)[segDfSubsetChr$end]))
                    irCNVR <- c(irCNVR,IRanges(start(inputChr)[start(cnvrRChr)],
                                                end(inputChr)[end(cnvrRChr)]))
                }
            }
        } else if (inputType=="DataMatrix"){
            ir <- IRanges(start=segDfSubset$start,end=segDfSubset$end)
            irCNVR <- IRanges(start=start(cnvrR),end=end(cnvrR))
        }

        rd <- GRanges(seqnames=segDfSubset$chr, ir,
                        seqinfo=seqinfo(grAllRegions), "sampleName"=sampleNames,
                        "median"=segDfSubset$median,"mean"=segDfSubset$mean,
                        "CN"=segDfSubset$CN)
        rd <- GenomeInfoDb::sortSeqlevels(rd)

        cnvr <- GRanges(seqnames=seqnames(cnvrR), irCNVR,
                        seqinfo=seqinfo(grAllRegions))
        cnvr <- GenomeInfoDb::sortSeqlevels(cnvr)
        values(cnvr) <- cnvrCN

        if (norm==2){
            r@normalizedData <- X.viz
        } else {
            r@normalizedData <- X.norm
        }
        r@localAssessments <- sINI
        r@individualCall   <- callsS
        r@iniCall          <- INI
        r@cnvs             <- rd
        r@cnvr             <- cnvr

        if (inputType=="GRanges"){
            irS <- IRanges()
            for (chrom in chrOrder){
                chrIdx <- chrDf[chrom,1]:chrDf[chrom,2]
                inputChr <- input[chrIdx]
                segDfChr <- subset(segDf,chr==chrom)
                if (nrow(segDfChr) >0 ){
                    irS <- c(irS, IRanges(start(inputChr)[segDfChr$start],
                                            end(inputChr)[segDfChr$end]))
                }
            }
            r@segmentation <- GRanges(seqnames=segDf$chr,
                                        irS, seqinfo=seqinfo(grAllRegions),
                                        "sampleName"=segDf$sample,
                                        "median"=segDf$median,
                                        "mean"=segDf$mean,"CN"=segDf$CN)
            r@segmentation <- GenomeInfoDb::sortSeqlevels(r@segmentation)

        } else if (inputType=="DataMatrix"){
            r@segmentation <- GRanges(seqnames=segDf$chr,
                                        IRanges(segDf$start,segDf$end),
                                        seqinfo=seqinfo(grAllRegions),
                                        "sampleName"=segDf$sample,
                                        "median"=segDf$median,
                                        "mean"=segDf$mean,"CN"=segDf$CN)
            r@segmentation <- GenomeInfoDb::sortSeqlevels(r@segmentation)
        }

        params <- append(params, list(badROI))
        names(params)[length(params)] <- "badROI"
        
        r@gr                <- grAllRegions
        r@posteriorProbs    <- post
        r@params            <- params
        r@integerCopyNumber <- CN
        r@sampleNames       <- colnames(X)
        message("Finished genepanelcn.mops")
        return(r)
    } else {
        message(paste("No CNVs detected. Try changing \"normalization\",",
                                    "\"priorImpact\" or \"I\"."))
        # Assembly of result object
        r <- new("CNVDetectionResult")
        return(r)
    }
}
