## Install the following packages: 
# CRAN packages: xlsx
# Follow the instructions at the Bioconductor website (https://bioconductor.org/)
# Search and install Bioconductor packages: GEOquery

## Load packages
library(Biobase)
library(GEOquery)
#library(ArrayExpress)
#library(xlsx)

## Set working directory
setwd("GOS_DataSets")
dir.create("GOS_DataSets/GEO") # here are store the dataset files downloaded by GEOquery package
dir.create("GOS_DataSets/Manual") # here are stored the dataset files downloaded manually from Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/)

## Download published preeclampsia (PE) datasets 
peNames <- c("GSE30186", "GSE35574", "GSE44711", "GSE60438", "GSE6573", "GSE73374", "GSE94643")
peList <- list()
for (i in 1:length(peNames)){
  peList[[i]] <- getGEO(peNames[i], GSEMatrix=TRUE, destdir = "GOS_DataSets/GEO", AnnotGPL=TRUE)
}
names(peList) <- peNames
## Check expression matrices
lapply(peList, function(x) exprs(x[[1]])[1:5,1:5])
## Get the methods of preprocessing
unlist(lapply(peList, function(x) as.character(x[[1]]@phenoData@data$data_processing[1])))
## Get probe annotation
peAnno <- lapply(peList, function(x) fData(x[[1]]))
lapply(peAnno, function(x) x[1:5,1:5])
## Get sample metadata
pePhen <- lapply(peList, function(x) pData(x[[1]]))
lapply(pePhen, function(x) x[1:5,1:5])
#for(i in 1:length(peList)){
#  write.xlsx(pePhen[[i]], file="PE_GSE_SampleInfo.xlsx", sheetName=names(pePhen)[i], 
#             col.names=TRUE, row.names=TRUE, append=TRUE)
#}

## Download published preterm birth (PB) and intrauterine growth restriction (IUGR) datasets
pbiuNames <- c("GSE35574", "GSE24129", "GSE73685")
pbiuList <- list()
for (i in 1:length(pbiuNames)){
  pbiuList[[i]] <- getGEO(pbiuNames[i], GSEMatrix=TRUE, destdir = "GOS_DataSets/GEO", AnnotGPL=TRUE)
}
names(pbiuList) <- pbiuNames
## Check expression matrices
lapply(pbiuList, function(x) exprs(x[[1]])[1:5,1:5])
## Get the methods of preprocessing
unlist(lapply(pbiuList, function(x) as.character(x[[1]]@phenoData@data$data_processing[1])))
## Get probe annotation
pbiuAnno <- lapply(pbiuList, function(x) fData(x[[1]]))
lapply(pbiuAnno, function(x) x[1:5,1:5])
## Get sample metadata
pbiuPhen <- lapply(pbiuList, function(x) pData(x[[1]]))
lapply(pbiuPhen, function(x) x[1:5,1:5])
#for(i in 1:length(pbiuList)){
#  write.xlsx(pbiuPhen[[i]], file="PB_IUGR_GSE_SampleInfo.xlsx", sheetName=names(pbiuPhen)[i], 
#             col.names=TRUE, row.names=TRUE, append=TRUE)
#}

## Save objects in R data file
save.image("Meta_GEOquery.RData")
