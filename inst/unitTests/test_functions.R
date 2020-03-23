test.getWindows <- function() {
    bed <- list.files(system.file("extdata", package = "panelcn.mops"),
                      pattern = ".bed$", full.names = TRUE)
    newCountWindows <- getWindows(bed)
    
    data(panelcn.mops)
    checkEquals(newCountWindows, countWindows,
                checkNames = FALSE)
}

test.runPanelcnMops <- function() {
    data(panelcn.mops)
    
    selectedGenes <- "ATM"
    
    XandCB <- test
    elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
                                     elementMetadata(control))
    newResultlist <- runPanelcnMops(XandCB, countWindows = countWindows, 
                                 selectedGenes = selectedGenes)
    checkEquals(integerCopyNumber(newResultlist[[1]]), 
                integerCopyNumber(resultlist[[1]]),
                checkNames = FALSE)
}

test.createResultTable <- function() {
    data(panelcn.mops)
    
    selectedGenes <- "ATM"
    
    XandCB <- test
    elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), 
                                     elementMetadata(control))
    sampleNames <- colnames(elementMetadata(test))
    newResulttable <- createResultTable(resultlist = resultlist, 
                                        XandCB = XandCB, 
                                        countWindows = countWindows, 
                                        selectedGenes = selectedGenes, 
                                        sampleNames = sampleNames)
    
    cns <- as.factor(c(rep("CN2", 60), rep("CN3", 2)))
    
    checkEquals(as.factor(newResulttable[[1]]$CN), cns,
                checkNames = FALSE)
}

