# panelcn.MOPS - CNV detection tool for targeted NGS panel data

The *panelcn.mops* R package is based on the *cn.mops* package and allows to detect copy number variations (CNVs) from targeted NGS panel data. Please visit http://www.bioinf.jku.at/software/panelcnmops/index.html for additional information.

## Installation:
1. install *cn.mops* from bioconductor.org:
   ```
   ## try http:// if https:// URLs are not supported  
   source("https://bioconductor.org/biocLite.R")  
   biocLite("cn.mops")  
   ```
2. install *panelcn.mops* from github e.g. like this:
   ```
   install.packages("devtools")  
   devtools::install_github("bioinf-jku/panelcn.mops")  
   ```