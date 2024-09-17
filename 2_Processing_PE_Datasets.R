## Install the following packages:
# - Bioconductor packages: GEOquery, affy, oligo, beadarray, limma, lumi, biomaRt
# - Bioconductor packages for microarray annotation: illuminaHumanv2.db, illuminaHumanv3.db, illuminaHumanv4.db, hgu133plus2.db, pd.hugene.2.0.st, pd.hugene.1.0.st.v1

## Load packages
library(GEOquery)
library(affy)
library(oligo)
library(beadarray)
library(limma)
library(lumi)
library(sva)
library(xlsx)
library(biomaRt)

## Load annotation packages
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

## Preeclampsia (PE) datasets

peFin <- list()

## GSE30186
tmp = tmp0 <- as(peList$GSE30186, "ExpressionSetIllumina")
tmp <- readBeadSummaryData(dataFile = "Manual/GSE30186_non_normalized.txt",
                           ProbeID = "ID_REF", skip = 0, qc.skip = 0,
                           illuminaAnnotation = "illuminaHumanv4.db")
fData(tmp) <- fData(tmp0)[featureNames(tmp), ]
View(pData(tmp)); View(pData(tmp0))
sampleNames(tmp) <- sampleNames(tmp0)
pData(tmp) <- data.frame(pData(tmp), pData(tmp0))
exprs(tmp) <- apply(exprs(tmp), 2, function(x) x + abs(min(x)) + 1)
tmp2 <- normaliseIllumina(tmp, method = "quantile", transform = "log2")
summary(exprs(tmp2))
plotSampleRelation(exprs(tmp2), method = "mds") # looks fine
peFin$GSE30186 <- tmp2
rm(tmp, tmp2)

## GSE35574
tmp <- as(peList$GSE35574, "ExpressionSetIllumina")
par(mfrow = c(1,2))
plotSampleRelation(exprs(tmp), method = "mds", col = c("red","blue")[as.numeric(factor(pData(tmp)$characteristics_ch1.5))]) # here we can see a clear batch effect, must be treated as two independent datasets
plotSampleRelation(exprs(tmp), method = "mds", col = c("red","blue")[as.numeric(factor(pData(tmp)$characteristics_ch1.1))])
dev.off()
table(pData(tmp)$characteristics_ch1.5, pData(tmp)$characteristics_ch1.1)
sampleNames(tmp)
# Get rid of duplicated samples
x <- c("GSM871048", "GSM871058", "GSM871064", "GSM871067")
table(sampleNames(tmp) %in% x)
tmp <- tmp[, !(sampleNames(tmp) %in% x)]
table(pData(tmp)$'classification:ch1')
tmp <- tmp[, pData(tmp)$'classification:ch1' != "IUGR"]
tmp2 <- tmp[, pData(tmp)$'characteristics_ch1.5' == "batch: A"]
table(as.character(pData(tmp)$'classification:ch1')); table(as.character(pData(tmp2)$'classification:ch1'))
peFin$GSE35574_A <- tmp2
tmp2 <- tmp[, pData(tmp)$'characteristics_ch1.5' == "batch: B"]
table(as.character(pData(tmp)$'classification:ch1')); table(as.character(pData(tmp2)$'classification:ch1'))
peFin$GSE35574_B <- tmp2
rm(tmp, tmp2)

## GSE44711
tmp = tmp0 <- as(peList$GSE44711, "ExpressionSetIllumina")
tmp <- readBeadSummaryData(dataFile = "Manual/GSE44711_non-normalized_data.txt",
                           ProbeID = "PROBE_ID", skip = 0, qc.skip = 0)
fData(tmp) <- fData(tmp0)[featureNames(tmp),]; fData(tmp)[1:5,1:8]
View(pData(tmp)); View(pData(tmp0))
sampleNames(tmp) <- sampleNames(tmp0)
pData(tmp) <- data.frame(pData(tmp), pData(tmp0))
tmp2 <- normaliseIllumina(tmp, method = "quantile", transform = "log2")
peFin$GSE44711 <- tmp2
rm(tmp, tmp0, tmp2)

## GSE60438
tmp <- as(peList$GSE60438, "ExpressionSetIllumina")
peFin$GSE60438 <- tmp
rm(tmp)

## GSE6573
tmp = tmp0 <- peList$GSE6573
x <- unlist(lapply(strsplit(as.character(pData(tmp)$title), " "), function(x) x[2]))
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(x))+1)
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$source_name_ch1))+1)
tmpaffy <- ReadAffy(filenames = list.celfiles(path="Manual/GSE6573_RAW", full.names=TRUE))
tmp <- affy::rma(tmpaffy)
fData(tmp) <- fData(tmp0)[featureNames(tmp),]; fData(tmp)[1:5,1:8]
sampleNames(tmp)  <-  sampleNames(tmp0)
pData(tmp) <- pData(tmp0)
table(pData(tmp)$source_name_ch1)
tmp <- tmp[, pData(tmp)$source_name_ch1 != "fat tissue"]
dim(tmp)
peFin$GSE6573 <- tmp 
rm(tmp, tmp0, tmpaffy, x)

## GSE73374
tmp = tmp0 <- peList$GSE73374
tmpaffy <- read.celfiles(filenames = list.celfiles(path="Manual/PE/GSE73374_RAW", full.names=TRUE),
                         pkgname = "pd.hugene.2.0.st")
tmp <- oligo::rma(tmpaffy, target="core")
featureData(tmp) <- getNetAffx(tmp, "transcript")
fData(tmp) <- data.frame(fData(tmp),fData(tmp0)[featureNames(tmp),])
sampleNames(tmp)  <-  sampleNames(tmp0)
pData(tmp) <- pData(tmp0)
summary(exprs(tmp)[,1:5]); summary(exprs(tmp0)[,1:5])
par(mfrow=c(1,2))
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$'characteristics_ch1'))+1)
plotSampleRelation(exprs(tmp0), method = "mds", col = as.numeric(factor(pData(tmp)$'characteristics_ch1'))+1)
dev.off()
View(pData(tmp))
table(pData(tmp)$characteristics_ch1)
peFin$GSE73374 <- tmp 
rm(tmp, tmp0, tmpaffy)

## GSE94643
tmp = tmp0 <- peList$GSE94643
tmpaffy <- read.celfiles(filenames = list.celfiles(path = "Manual/GSE94643_RAW", full.names = TRUE),
                         pkgname = "pd.hugene.2.0.st")
tmp <- rma(tmpaffy, target = "core")
featureData(tmp) <- getNetAffx(tmp, "transcript")
fData(tmp) <- data.frame(fData(tmp),fData(tmp0)[featureNames(tmp),])
sampleNames(tmp)  <-  sampleNames(tmp0)
pData(tmp) <- pData(tmp0)
summary(exprs(tmp)[,1:5]); summary(exprs(tmp0)[,1:5])
par(mfrow=c(1,2))
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$'characteristics_ch1'))+1)
plotSampleRelation(exprs(tmp0), method = "mds", col = as.numeric(factor(pData(tmp)$'characteristics_ch1'))+1)
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$'characteristics_ch1'))+1)
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(tmp)$'source_name_ch1'))+1)
dev.off()
View(pData(tmp))
table(pData(tmp)$characteristics_ch1)
table(pData(tmp)$characteristics_ch1, pData(tmp)$source_name_ch1)
peFin$GSE94643 <- tmp 
rm(tmp, tmp0, tmpaffy)

## In-home dataset (Russian and Yakut women with or without preeclampsia)
library("illuminaHumanv4.db")
load("PE_NIIMG_Data.RData") # This data is available upon request (e-mail to )
annotation(BSData) <- "illuminaHumanv4.db"
fData(BSData)$PROBE_ID <- unlist(mget(as.character(fData(BSData)$ProbeID), revmap(illuminaHumanv4ARRAYADDRESS), ifnotfound = NA))
BSData <- BSData[!(is.na(as.character(fData(BSData)$PROBE_ID))),]
featureNames(BSData) <- as.character(fData(BSData)$PROBE_ID)
tmp <- mget(featureNames(BSData), illuminaHumanv4REFSEQ, ifnotfound = NA)
fData(BSData)$RefSeq_ID <- as.character(fData(peFin[[1]])[featureNames(BSData),'RefSeq_ID'])
fData(BSData)$Entrez_Gene_ID <- as.character(fData(peFin[[1]])[featureNames(BSData),'Entrez_Gene_ID'])
fData(BSData)$Definition <- as.character(fData(peFin[[1]])[featureNames(BSData),'Definition'])
tmp <- BSData
plotSampleRelation(exprs(tmp), method = "mds", col = as.numeric(factor(pData(BSData)$Nationality))+1) # batch by ethnicity
head(pData(BSData))
# Keep Russian samples apart
tmp <- BSData[, pData(BSData)$Nationality == "R"]
tmp <- normaliseIllumina(tmp, method = "quantile", transform = "log2")
summary(exprs(tmp))
pData(tmp)
table(pData(tmp)$Diagnosis)
peFin$OUR_R <- tmp
# Keep Yakut samples apart
tmp <- BSData[, pData(BSData)$Nationality == "Y"]
tmp <- normaliseIllumina(tmp, method = "quantile", transform = "log2")
summary(exprs(tmp))
pData(tmp)
table(pData(tmp)$Diagnosis)
peFin$OUR_Y <- tmp

## Get PE sample information
pePhen <- lapply(peFin, function(x) pData(x))
lapply(pePhen, function(x) x[1:5,1:5])
for(i in 1:length(pePhen)){
  write.xlsx(pePhen[[i]], file="PE_Samples.xlsx", sheetName=names(pePhen)[i], 
             col.names=TRUE, row.names=TRUE, append=TRUE)
}

## Get PE probe annotation
lapply(peFin, function(x) dim(fData(x)))
peAnno <- lapply(peFin, function(x) fData(x))
for(i in 1:length(peAnno)){
  write.xlsx(peAnno[[i]][10000:10100,], file="PE_Annotation.xlsx", sheetName=names(peAnno)[i], 
             col.names=TRUE, row.names=TRUE, append=TRUE)
}

## Load dataset information
gseInfo <- read.delim("Datasets_Info.tsv", na.strings = "")
## Choose only PE datasets
peInfo <- gseInfo[gseInfo$Group == "PE",]

## Save results
save(peFin, peAnno, pePhen, peInfo, file = "Meta_PE_Fin.RData")

## Manual merging of PE datasets

## Load sample information for PE dataset
peSampleInfo <- read.delim("PE_Sample_Info.tsv")
rownames(peSampleInfo) <- peSampleInfo$Sample

## Calculate MAD for each probe in each dataset
peMad <- lapply(peFin, function(eset) return(apply(exprs(eset), 1, mad)))
lapply(peMad, head)
unlist(lapply(peMad, length))
unlist(lapply(peFin, function(x) dim(x)[1]))

## Save results
save(peFin, peAnno, pePhen, peInfo, file = "Meta_PE_Fin.RData")

## Get access to BioMart database
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
martFilters <- listFilters(mart)
attributes <-  listAttributes(mart)

## Annotate datasets with Ensembl IDs and HGNC Symbols using BioMart
peMartList <- list()
for(i in 1:nrow(peInfo)){
  myfilter <- as.character(gseInfo$Platform_Mart)[i]
  values <- unique(as.character(rownames(peAnno[[i]])))
  peMartList[[i]] <- getBM(attributes = c(myfilter, 'ensembl_gene_id', 'hgnc_symbol'), 
               filters = myfilter, 
               values = values, 
               mart = mart)
  #rownames(selMart) <- as.character(selMart[,1])
}; names(peMartList) <- names(peAnno); rm(i, myfilter, values)
lapply(peMartList, head)
lapply(peMartList, nrow)
lapply(peMartList, function(x) table(duplicated(x[,1])))

## Add Ensembl IDs to annotations 
for(i in 1:length(peAnno)){
  x <- peMartList[[i]]
  x <- x[!duplicated(x[,1]), ]
  peAnno[[i]]$BM_Symbol = peAnno[[i]]$BM_Ensembl_ID <- rep(NA, nrow(peAnno[[i]]))
  peAnno[[i]][as.character(x[,1]), 'BM_Ensembl_ID'] <- as.character(x[,2])
  peAnno[[i]][as.character(x[,1]), 'BM_Symbol'] <- as.character(x[,3])
}; rm(i,x)

## Get shared Ensembl IDs
peShared <- Reduce(intersect, lapply(peMartList, function(x) return(as.character(x[,'ensembl_gene_id']))))
length(peShared) ## 21598 IDs
peShared <- Reduce(intersect, lapply(peAnno, function(x) return(as.character(na.omit(x[,'BM_Ensembl_ID'])))))
length(peShared) ## 16781 is fine too
peSharedTab <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                     filters = 'ensembl_gene_id', values = peShared, mart = mart)
peSharedTab <- peSharedTab[!duplicated(peSharedTab$ensembl_gene_id),]
dim(peSharedTab)
rownames(peSharedTab) <- peSharedTab$ensembl_gene_id

### Create the list of expressionSets gathering only shared probes 
peFinShared <- list()
for(i in 1:length(peFin)){
  eset <- peFin[[i]]
  fData(eset) <- peAnno[[i]]
  eset <- eset[order(peMad[[i]], decreasing = T), ]
  refseq <- fData(eset)$BM_Ensembl_ID
  sel <- which((refseq %in% peShared) & (!duplicated(refseq)))
  eset <- eset[sel,]
  featureNames(eset) <- as.character(fData(eset)$BM_Ensembl_ID)
  eset <- eset[peShared,]
  fData(eset) <- peSharedTab
  pData(eset) <- peSampleInfo[sampleNames(eset),]
  peFinShared[[i]] <- eset
}
rm(eset, refseq, sel)
names(peFinShared) <- names(peFin)
lapply(peFinShared, dim)
lapply(peFinShared, function(x) featureNames(x)[1:5])
lapply(peFinShared, function(x) head(fData(x)))

### Create and save expression matrix, sample information and probe annotation
pe.exprs <- do.call(cbind, lapply(peFinShared, exprs))
pe.pd <- do.call(rbind, lapply(peFinShared, pData))
rownames(pe.pd) <- as.character(pe.pd$Sample)
pe.pd$Group <- factor(pe.pd$Group)
pe.pd$GroupShort <- pe.pd$Group
levels(pe.pd$GroupShort) <- c("CTRL", "PE")
pe.pd$Dataset <- factor(pe.pd$Dataset)
pe.fd <- fData(peFinShared[[1]])
save(pe.exprs, pe.pd, pe.fd, file = "Meta_PE_BioMart.RData")
