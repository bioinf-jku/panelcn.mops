## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(

)

## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----echo=FALSE, message=FALSE---------------------------------
options(width=65)
set.seed(0)
library(panelcn.mops)
panelcn.mopsVersion <- packageDescription("panelcn.mops")$Version

## ----echo=TRUE-------------------------------------------------
library(panelcn.mops)

## ----eval=FALSE------------------------------------------------
#  bed <- "Genes_part.bed"
#  countWindows <- getWindows(bed)

## ----eval=FALSE------------------------------------------------
#  testbam <- "SAMPLE1.bam"
#  test <- countBamListInGRanges(countWindows = countWindows,
#                                  bam.files = testbam, read.width = 150)

## ----eval=FALSE------------------------------------------------
#  
#  data(panelcn.mops)
#  
#  selectedGenes <- c("ATM")
#  
#  XandCB <- test
#  elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
#                                  elementMetadata(control))
#  
#  resultlist <- runPanelcnMops(XandCB,
#                              testiv = 1:ncol(elementMetadata(test)),
#                              countWindows = countWindows,
#                              selectedGenes = selectedGenes)

## ----echo=FALSE,results='hide'---------------------------------
data(panelcn.mops)

XandCB <- test
elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
                                elementMetadata(control))
selectedGenes <- "ATM"

## --------------------------------------------------------------
sampleNames <- colnames(elementMetadata(test))
resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, 
                                    countWindows = countWindows, 
                                    selectedGenes = selectedGenes, 
                                    sampleNames = sampleNames)

(tail(resulttable[[1]]))

## ----eval=FALSE------------------------------------------------
#  plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1],
#              countWindows = countWindows,
#              selectedGenes = selectedGenes, showGene = 1)

## ----fig.keep='none',echo=FALSE,results='hide'-----------------
data(panelcn.mops)
sampleNames <- colnames(elementMetadata(test))
selectedGenes <- "ATM"

pdf("001.pdf", width = 10)
plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1], 
            countWindows = countWindows, selectedGenes = selectedGenes, 
            showGene = 1)
dev.off()

## --------------------------------------------------------------
bed <- list.files(system.file("extdata", package = "panelcn.mops"),
                        pattern = ".bed$", full.names = TRUE)
countWindows <- getWindows(bed)

## ----echo=FALSE------------------------------------------------
bed <- list.files(system.file("extdata", package = "panelcn.mops"),
                        pattern = ".bed$", full.names = TRUE)
write.table(head(read.table(bed)), row.names = FALSE, col.names = FALSE, 
            quote = FALSE)

## ----eval=FALSE------------------------------------------------
#  testbam <- "SAMPLE1.bam"
#  test <- countBamListInGRanges(countWindows = countWindows,
#                                  bam.files = testbam, read.width = 150)

## --------------------------------------------------------------
(test)

## --------------------------------------------------------------
splitROIs(bed, "newBed.bed")

## ----eval=FALSE------------------------------------------------
#  
#  data(panelcn.mops)
#  
#  XandCB <- test
#  elementMetadata(XandCB) <- cbind(elementMetadata(XandCB),
#                                  elementMetadata(control))
#  resultlist <- runPanelcnMops(XandCB, countWindows = countWindows)

## ----eval=FALSE------------------------------------------------
#  (str(resultlist[[1]]))

## ----eval=FALSE------------------------------------------------
#  help(CNVDetectionResult)

## --------------------------------------------------------------
integerCopyNumber(resultlist[[1]])[1:5]

## --------------------------------------------------------------
sampleNames <- colnames(elementMetadata(test))
resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, 
                                    countWindows = countWindows, 
                                    selectedGenes = selectedGenes, 
                                    sampleNames = sampleNames)

(tail(resulttable[[1]]))

## ----eval=FALSE------------------------------------------------
#  plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1],
#              countWindows = countWindows,
#              selectedGenes = selectedGenes, showGene = 1)

## ----eval=FALSE------------------------------------------------
#  toBibtex(citation("panelcn.mops"))

