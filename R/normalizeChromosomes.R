# Copyright (C) 2011 Klambauer Guenter
# <klambauer@bioinf.jku.at>

.statmod <- function(x,na.rm=FALSE) {
    if (na.rm){
        z <- table(as.vector(x[!is.na(x)]))
        r <- names(z)[z == max(z)]
        return(as.numeric(r)[1])
    } else {
        if (any(is.na(x))){return(NA)
        } else {
            z <- table(as.vector(x))
            r <- names(z)[z == max(z)]
            return(as.numeric(r)[1])
        }
    }
}


#' @title Normalization of NGS data.
#'
#' @description Normalize quantitative NGS data in order to make counts
#' comparable over samples, i.e., correcting for different library sizes or
#' coverages. Scales each samples' reads such that the coverage is even for
#' all samples after normalization.
#' @param X Matrix of positive real values, where columns are interpreted as
#' samples and rows as genomic regions. An entry is the read count of a sample
#' in the genomic region. Alternatively this can be a GRanges object
#' containing the read counts as values.
#' @param chr Character vector that has as many elements as "X" has rows. The
#' vector assigns each genomic segment to a reference sequence (chromosome).
#' @param normType Type of the normalization technique. Each samples'
#' read counts are scaled such that the total number of reads are comparable
#' across samples.
#' If this parameter is set to the value "mode",
#' the read counts are scaled such that each samples'
#' most frequent value (the "mode") is equal after normalization.
#' Accordingly for the other options are "mean","median","poisson", "quant",
#' and "mode".
#' Default = "poisson".
#' @param sizeFactor  By this parameter one can decide to how the size factors
#' are calculated.
#' Possible choices are the mean, median or mode coverage
#' ("mean", "median", "mode") or any quantile ("quant").
#' @param qu Quantile of the normType if normType is set to "quant".
#' Real value between 0 and 1. Default = 0.25.
#' @param quSizeFactor Quantile of the sizeFactor if sizeFactor is set to
#' "quant". 0.75 corresponds to "upper quartile normalization".
#' Real value between 0 and 1. Default = 0.75.
#' @param ploidy An integer value for each sample or each column in the read
#' count matrix. At least two samples must have a ploidy of 2.
#' Default = "missing".
#' @examples
#' data(cn.mops)
#' X.norm <- normalizeChromosomes(X)
#' @return A data matrix of normalized read counts with the same dimensions
#' as the input matrix X.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at}
#' @export


normalizeChromosomes <- function(X, chr, normType="poisson", sizeFactor="mean",
                                    qu=0.25, quSizeFactor=0.75, ploidy){
    if (!(normType %in% c("mean", "median","poisson", "mode", "quant"))){
        stop(paste("Set \"normType\" of normalization to \"mean\"",
                        "\"median\", \"quant\", \"mode\" or \"poisson\"."))
    }
    if (!(sizeFactor %in% c("mean","median","quant","mode"))){
        stop(paste("Set \"sizeFactor\" of normalization to \"mean\"",
                        "\"median\", \"quant\" or \"mode\"."))
    }
    input <- X
    if (!(is.matrix(input)|class(input)=="GRanges")){
        stop("Input data must be matrix or GRanges object!")
    }
    returnGRanges <- FALSE
    if(class(X)=="GRanges"){
        returnGRanges <- TRUE
        X <- IRanges::as.matrix(IRanges::values(X))
    }
    if (is.vector(X)){X <- matrix(X,nrow=1)}

    if (missing(chr) & is.matrix(input)){
        chr <- rep("undef",nrow(X))
    }
    if (missing(chr) & class(input)=="GRanges"){
        chr <- as.character(seqnames(input))
    }
    if (missing(ploidy)){
        ploidy <- rep(2,ncol(X))
    }
    if (any(ploidy!=as.integer(ploidy))){
        stop("Ploidy values must be integers!")
    }
    if (length(ploidy)!=ncol(X)){
        stop("Length of the ploidy vector does not match the number of",
                "columns of the read count matrix!")
    }
    ploidy <- as.integer(ploidy)
    if (length(unique(ploidy))==1) ploidy <- rep(2, ncol(X))
    if (!length(which(ploidy>=2))){
        stop("At least two diploid samples must be contained in the data.")
    }
    ploidy[ploidy==0] <- 0.05
    Xorig <- X

    # Sequencing data matrix
    # vector of chromosome - length equal to rows of X
    if (length(chr)!=nrow(X)){
        stop("Length of \"chr\" must be equal to number of rows of \"X\".")}
    chr <- (as.character(chr))

    YY <- matrix(0,nrow=nrow(Xorig),ncol=ncol(Xorig))

    ploidy2flag <- FALSE
    ploidy2median <- c()
    for (pp in unique(c(2,ploidy))){
        X <- Xorig[,ploidy==pp,drop=FALSE]
        globalsizeFactors <- colSums(X)

        if (any(is.na(globalsizeFactors))) stop("NAs in readcount matrix.")
        if (any(globalsizeFactors==0)) stop("Zero columns in readcount matrix.")

        if (ncol(X)==1){
            Y <- X
        } else {

            Y <- matrix(0,nrow=nrow(X),ncol=ncol(X))
            for (l in (unique(chr))){
                chrIdx <- which(chr==l)
                Ytmp <- X[chrIdx, ,drop=FALSE]

                if (all(Ytmp==0)){
                    Y[chrIdx, ] <- Ytmp
                } else {
                    idxSG <- apply(Ytmp,1,function(x) all(x<1))
                    Ytmp[idxSG, ] <- NA

                    if ((nrow(Ytmp)-length(which(idxSG))) > 1){

                        if (sizeFactor=="mean"){
                            sizeFactors <- colMeans(Ytmp,na.rm=TRUE)
                        } else if (sizeFactor=="median"){
                            sizeFactors <- apply(Ytmp, 2, median, na.rm=TRUE)
                        } else if (sizeFactor=="quant"){
                            sizeFactors <- apply(Ytmp, 2, quantile,
                                                probs=quSizeFactor,na.rm=TRUE)
                        } else if (sizeFactor=="mode"){
                            sizeFactors <- apply(Ytmp, 2, function(x)
                                                .statmod(x[x!=0],na.rm=TRUE))
                        }

                        if (any(is.na(sizeFactors))){
                            warning(paste("Normalization failed for
                                            reference sequence ", l,
                                            ". Using global sizeFactors!"))
                            sizeFactors <- globalsizeFactors
                        }

                        if (any(sizeFactors==0)){
                            stop(paste("Some normalization factors are zero!",
                                "Remove samples or chromosomes for which the",
                                "average read count is zero,",
                                "e.g. chromosome Y."))
                        }


                        if (normType=="mean"){
                            correctwiththis <-
                                mean(sizeFactors,na.rm=TRUE)/sizeFactors
                        } else if (normType=="median"){
                            correctwiththis <-
                                median(sizeFactors,na.rm=TRUE)/sizeFactors
                        } else if (normType=="quant"){
                            correctwiththis <- quantile(sizeFactors, probs=qu,
                                                        na.rm=TRUE)/sizeFactors
                        } else if (normType=="mode"){
                            asm <- .statmod(Ytmp[Ytmp!=0],na.rm=TRUE)
                            correctwiththis <- asm/sizeFactors
                        } else if (normType=="poisson"){
                            correctwiththis <-
                                mean(sizeFactors,na.rm=TRUE)/sizeFactors
                            YYtmp <- t(t(Ytmp)*correctwiththis)
                            v2m <- apply(YYtmp,1,var)/rowMeans(YYtmp)
                            uut <- quantile(v2m,probs=0.95,na.rm=TRUE)
                            dd <- density(v2m,na.rm=TRUE,from=0,
                                    to=uut,n=uut*100)
                            #mv2m <- median(v2m,na.rm=TRUE)
                            mv2m <- dd$x[which.max(dd$y)]
                            if (is.finite(mv2m)) {
                                correctwiththis <- correctwiththis*1/mv2m
                            }
                        }
                        #browser()

                    } else {
                        warning(paste("Normalization for reference sequence ",l,
                        "not applicable, because of low number of segments"))
                        correctwiththis <- rep(1,ncol(X))
                    }

                    if (any(!is.finite(correctwiththis))){
                        warning(paste("Normalization for reference sequence ",l,
                                        "not applicable, because at least one",
                                        "sample has zero reads."))
                        correctwiththis <- rep(1,ncol(X))
                    }


                    Ytmp <- t(t(Ytmp)*correctwiththis)
                    Ytmp[idxSG, ] <- 0

                    Y[chrIdx, ] <- Ytmp
                }
            }
        } # over chr

        if (!ploidy2flag){
            ploidy2flag <- TRUE
            ploidy2median <- median(Y[!idxSG, ],na.rm=TRUE)
        }

        if (pp!=2){
            mm <- median(Y[!idxSG, ],na.rm=TRUE)
            if (ploidy2median==0 & mm==0){
                YY[,ploidy==pp] <- Y*pp/2
            } else {
                YY[,ploidy==pp] <- Y*ploidy2median/mm*pp/2
            }
        } else{
            YY[,ploidy==pp] <- Y
        }
    }

    rownames(YY) <- rownames(Xorig)
    colnames(YY) <- colnames(Xorig)

    if (returnGRanges){
        values(input) <- YY
        return(input)
    } else {
        return(YY)
    }
}

