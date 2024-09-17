#### Load packages
library(GEOquery)
library(affy)
library(oligo)
library(beadarray)
library(limma)
library(lumi)
library(sva)
library(xlsx)
library(biomaRt)

## Annotation packages
x <- c("illuminaHumanv2.db", "illuminaHumanv3.db", "illuminaHumanv4.db", 
       "hgu133plus2.db", "pd.hugene.2.0.st", "pd.hugene.1.0.st.v1")
lapply(x, require, character.only = TRUE)

## Set working directory
setwd("GOS_DataSets")

## Load the downloaded datasets
load("Meta_GEOquery.RData")

## Tidy up and normalise the dataset 
# We use either pre-processed data downloaded by GEOquery or manually downloaded raw data. 
# We process the data up to log2-scaled and quantile normalized expression values.

## Preterm birth (PB) and intrauterine growth restriction (IUGR) datasets

pbiuFin <- list()

## GSE35574
tmp <- as(pbiuList$GSE35574, "ExpressionSetIllumina")
par(mfrow=c(1,2))
plotSampleRelation(exprs(tmp), method = "mds", col = c("red","blue")[as.numeric(factor(pData(tmp)$characteristics_ch1.5))]) # batch effect
plotSampleRelation(exprs(tmp), method = "mds", col = c("red","blue")[as.numeric(factor(pData(tmp)$characteristics_ch1.1))])
dev.off()
table(pData(tmp)$characteristics_ch1.5, pData(tmp)$characteristics_ch1.1)
sampleNames(tmp)
# Get rid of duplicated samples
x <- c("GSM871048", "GSM871058", "GSM871064", "GSM871067")
table(sampleNames(tmp) %in% x)
tmp <- tmp[, !(sampleNames(tmp) %in% x)]
table(pData(tmp)$'classification:ch1')
tmp <- tmp[, pData(tmp)$'classification:ch1' != "PE"]
tmp2 <- tmp[, pData(tmp)$'characteristics_ch1.5' == "batch: A"]
table(as.character(pData(tmp)$'classification:ch1')); table(as.character(pData(tmp2)$'classification:ch1'))
pbiuFin$GSE35574_A <- tmp2
tmp2 <- tmp[, pData(tmp)$'characteristics_ch1.5' == "batch: B"]
table(as.character(pData(tmp)$'classification:ch1')); table(as.character(pData(tmp2)$'classification:ch1'))
pbiuFin$GSE35574_B <- tmp2
rm(tmp, tmp2)

## GSE24129
tmp = tmp0 <- pbiuList$GSE24129
tmpaffy <- read.celfiles(filenames = list.celfiles(path = "Manual/GSE24129_RAW", full.names = TRUE), pkgname = "pd.hugene.1.0.st.v1")
tmp <- rma(tmpaffy, target="core")
featureData(tmp) <- getNetAffx(tmp, "transcript")
fData(tmp) <- data.frame(fData(tmp),fData(tmp0)[featureNames(tmp),])
sampleNames(tmp)  <-  sampleNames(tmp0)
pData(tmp) <- pData(tmp0)
sampleNames(tmp) <- sampleNames(tmp0)
summary(exprs(tmp)[,1:5]); summary(exprs(tmp0)[,1:5])
par(mfrow=c(1,2))
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$'disease'))+1)
plotSampleRelation(exprs(tmp0), method = "mds", col = as.numeric(factor(pData(tmp0)$'disease'))+1)
plotDensities(tmp, legend = FALSE); plotDensities(tmp0, legend = FALSE) ## Rather good normalization from authors
dev.off()
table(pData(tmp)$"disease")
tmp <- tmp[, pData(tmp)$"disease" != "Pre-eclampsia"]
pbiuFin$GSE24129 <- tmp 
rm(tmp, tmp0, tmpaffy)

## GSE73685 
tmp = tmp0 <- pbiuList$GSE73685
boxplot(tmp)
plotDensities(tmp, legend = FALSE)
tmpaffy <- read.celfiles(filenames = list.celfiles(path = "Manual/GSE73685_RAW", full.names = TRUE), pkgname = "pd.hugene.1.0.st.v1")
tmp <- rma(tmpaffy, target = "core")
featureData(tmp) <- getNetAffx(tmp, "transcript")
fData(tmp) <- data.frame(fData(tmp),fData(tmp0)[featureNames(tmp),])
sampleNames(tmp)  <-  sampleNames(tmp0)
pData(tmp) <- pData(tmp0)
summary(exprs(tmp)[,1:5]); summary(exprs(tmp0)[,1:5])
par(mfrow=c(1,2))
gr <- factor(pData(tmp)$'tissue type:ch1')
table(gr)
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(gr)+1) # there is a grouping of samples
plotSampleRelation(exprs(tmp0), method = "mds", col = as.numeric(gr)+1) # the grouping is tissue-specific
dev.off()
# we keep only Chorion, Placenta and Decidua for further analysis
tmp2 <- tmp[, gr %in% c("Chorion", "Decidua", "Placenta")]
# And we need only preterm (PL, PNL) and normal term (TL,TNL) samples
table(pData(tmp2)$'outcome:ch1')
tmp2 <- tmp2[, !(pData(tmp2)$'outcome:ch1' %in% c("pPROM no labor", "pPROM with labor"))]
pbiuFin$GSE73685 <- tmp2
rm(tmp2, tmp, tmp0, tmpaffy, gr)

## Inspect the final list
lapply(pbiuFin, function(x) summary(exprs(x[, 1:3])))

## Get phenotype data
pbiuPhen <- lapply(pbiuFin, function(x) pData(x))
lapply(pbiuPhen, function(x) x[1:5,1:5])
for(i in 1:length(pbiuPhen)){
  write.xlsx(pbiuPhen[[i]], file = "PBIU_Samples.xlsx", sheetName = names(pbiuPhen)[i], 
             col.names = TRUE, row.names = TRUE, append = TRUE)
}

## Get annotation data
lapply(pbiuFin, function(x) dim(fData(x)))
pbiuAnno <- lapply(pbiuFin, function(x) fData(x))
lapply(pbiuAnno, function(x) x[10000:10010,1:8])
for(i in 1:length(pbiuAnno)){
  write.xlsx(pbiuAnno[[i]][10000:10100,], file = "PBIU_Annotation.xlsx", sheetName = names(pbiuAnno)[i], 
             col.names = TRUE, row.names = TRUE, append = TRUE)
}

gseInfo <- read.delim("Datasets_Info_2020.txt")
pbiuInfo <- gseInfo[gseInfo$Group != "PE",]

## Save results
save(pbiuFin, pbiuAnno, pbiuPhen, pbiuInfo, file = "Meta_PBIU_Fin.RData")

## Annotate datasets with Ensembl IDs and HGNC Symbols using BioMart
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
martFilters <- listFilters(mart)
attributes <-  listAttributes(mart)
pbiuMartList <- list()
for(i in 1:length(pbiuAnno)){
  myfilter <- as.character(pbiuInfo$Platform_Mart)[i]
  values <- unique(as.character(rownames(pbiuAnno[[i]])))
  pbiuMartList[[i]] <- getBM(attributes=c(myfilter, 'ensembl_gene_id', 'hgnc_symbol'), 
                           filters = myfilter, 
                           values = values, 
                           mart = mart)
}; names(pbiuMartList) <- names(pbiuAnno); rm(i, myfilter, values)
lapply(pbiuMartList, head)
lapply(pbiuMartList, dim)
lapply(pbiuMartList, function(x) table(duplicated(x[,1])))
## Add Ensembl IDs to annotations 
for(i in c(1:length(pbiuAnno))){
  x <- pbiuMartList[[i]]
  x <- x[!duplicated(x[,1]), ]
  pbiuAnno[[i]]$BM_Symbol = pbiuAnno[[i]]$BM_Ensembl_ID <- rep(NA, nrow(pbiuAnno[[i]]))
  pbiuAnno[[i]][as.character(x[,1]), 'BM_Ensembl_ID'] <- as.character(x[,2])
  pbiuAnno[[i]][as.character(x[,1]), 'BM_Symbol'] <- as.character(x[,3])
}; rm(i,x)
lapply(pbiuAnno, colnames)

## Calculate MAD for each probe in each dataset
pbiuMad <- lapply(pbiuFin, function(eset) return(apply(exprs(eset), 1, mad)))
lapply(pbiuMad, head)
lapply(pbiuMad, summary)
unlist(lapply(pbiuMad, length))

## Manual merging of IUGR datasets
## Load sample info for IUGR
iugrSampleInfo <- read.delim("IUGR_Sample_Info.tsv")
rownames(iugrSampleInfo) <- iugrSampleInfo$Sample
iugrSampleInfo$GroupShort <- factor(iugrSampleInfo$Group)
levels(iugrSampleInfo$GroupShort) <- c("CTRL", "IUGR")
## Get shared Ensembl IDs
iugrShared <- Reduce(intersect, lapply(pbiuAnno[which(pbiuInfo$Group == "IUGR")], 
                                       function(x) return(as.character(na.omit(x[,'BM_Ensembl_ID'])))))
length(iugrShared)
iugrSharedTab <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                     filters = 'ensembl_gene_id', values = iugrShared, mart = mart)
iugrSharedTab <- iugrSharedTab[!duplicated(iugrSharedTab$ensembl_gene_id),]
rownames(iugrSharedTab) <- iugrSharedTab$ensembl_gene_id
## Create the list of expressionSets gathering only shared probes 
iugrFinShared <- list()
for(i in 1:length(pbiuFin)){
  eset <- pbiuFin[[i]]
  fData(eset) <- pbiuAnno[[i]]
  eset <- eset[order(pbiuMad[[i]], decreasing=T), ]
  refseq <- fData(eset)$BM_Ensembl_ID
  sel <- which((refseq %in% iugrShared) & (!duplicated(refseq)))
  eset <- eset[sel,]
  featureNames(eset) <- as.character(fData(eset)$BM_Ensembl_ID)
  eset <- eset[iugrShared,]
  fData(eset) <- iugrSharedTab
  pData(eset) <- iugrSampleInfo[sampleNames(eset),]
  iugrFinShared[[i]] <- eset
}
rm(eset, refseq, sel)
names(iugrFinShared) <- names(pbiuFin)[1:3]
lapply(iugrFinShared, dim)
lapply(iugrFinShared, function(x) head(exprs(x)[, 1:3]))
lapply(iugrFinShared, function(x) head(pData(x)))
## Create expression matrix, phenoData and annotation
iugr.exprs <- do.call(cbind, lapply(iugrFinShared, exprs))
iugr.pd <- do.call(rbind, lapply(iugrFinShared, pData))
rownames(iugr.pd) <- as.character(iugr.pd$Sample)
iugr.fd <- fData(iugrFinShared[[1]])

## Prepare PB dataset
## Load sample info for PB datasets
pbSampleInfo <- read.delim("PB_Sample_Info.tsv")
rownames(pbSampleInfo) <- pbSampleInfo$Sample
## Get shared Ensembl IDs
pbShared <- Reduce(intersect, lapply(pbiuAnno[which(pbiuInfo$Group == "PB")], 
                                     function(x) return(as.character(na.omit(x[,'BM_Ensembl_ID'])))))
length(pbShared) 
pbSharedTab <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                     filters = 'ensembl_gene_id', values = pbShared, mart = mart)
pbSharedTab <- pbSharedTab[!duplicated(pbSharedTab$ensembl_gene_id),]
rownames(pbSharedTab) <- pbSharedTab$ensembl_gene_id
pbFinShared <- pbiuFin[["GSE73685"]]
fData(pbFinShared) <- pbiuAnno[["GSE73685"]]
eset <- pbFinShared[order(pbiuMad[["GSE73685"]], decreasing = T),]
refseq <- fData(pbFinShared)$BM_Ensembl_ID
sel <- which((refseq %in% pbShared) & (!duplicated(refseq)))
pbFinShared <- pbFinShared[sel,]
featureNames(pbFinShared) <- as.character(fData(eset)$BM_Ensembl_ID)
pbFinShared <- pbFinShared[pbShared,]
fData(pbFinShared) <- pbSharedTab
pData(pbFinShared) <- pbSampleInfo[sampleNames(pbFinShared),]
rm(refseq, sel)
dim(pbFinShared)
head(exprs(pbFinShared)[, 1:3])
head(pData(pbFinShared))
## Create  expression matrices, sample information and probe annotation
pb.exprs <- exprs(pbFinShared)
pb.pd <- pData(pbFinShared)
rownames(pb.pd) <- as.character(pb.pd$Sample)
pb.pd$Group <- factor(pb.pd$Group)
pb.pd$GroupShort <- pb.pd$Group
levels(pb.pd$GroupShort) <- c("CTRL", "PB")
pb.pd$Dataset <- factor(pb.pd$Dataset)
pb.fd <- fData(pbFinShared)

## Save results
save(iugr.exprs, iugr.pd, iugr.fd, pb.exprs, pb.pd, pb.fd, file = "Meta_IUGR_PB_BioMart.RData")
