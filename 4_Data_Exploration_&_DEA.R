## Load packages
library(limma)
library(sva)
library(xlsx)
library(RColorBrewer)

## Set color palettes
mypal1 <- brewer.pal(8, "Set1")[-6]
mypal2 <- brewer.pal(7, "Set2")
mypal3 <- c(mypal1, mypal2)

## Set working directory
setwd("GOS_DataSets")

## Load the processed datasets
load("Meta_PE_BioMart.RData")
load("Meta_PE_IUGR_PB_BioMart.RData")

#### Reanalysis of preeclampsia (PE) samples ####

table(pe.pd$Tissue, pe.pd$Group)
pe.pd$Dataset <- factor(pe.pd$Dataset)
pe.pd$Tissue <- factor(pe.pd$Tissue)
pe.pd$Country <- factor(pe.pd$Country)
## MDS plots
tmp <- pe.exprs
plotMDS(tmp, labels = pe.pd$Country, 
        col=mypal[2:1][pe.pd$Group])
plotMDS(tmp, top = nrow(tmp), labels = pe.pd$Dataset, 
        col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country])
legend('bottomleft', legend = levels(pe.pd$Country), pch=15, col = mypal[1:length(unique(pe.pd$Country))])

plotMDS(tmp, labels = pe.pd$Dataset,
        col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country])
plotMDS(tmp, top = nrow(tmp), labels = pe.pd$GroupShort, pch=19,
        col=mypal[1:length(unique(pe.pd$Tissue))][pe.pd$Tissue])
legend('bottomleft', legend = levels(pe.pd$Tissue), pch = 15, col = mypal[1:length(unique(pe.pd$Tissue))])
par(mfrow=c(1,2))
tmp.tsne <- Rtsne(tmp[, pe.pd$Group == "Control"], check_duplicates = FALSE, pca = TRUE, perplexity = 10, theta = 0.5, dims = 2)
plot(tmp.tsne$Y, main="Control")
tmp.tsne <- Rtsne(tmp[, pe.pd$Group == "Preeclampsia"], check_duplicates = FALSE, pca = TRUE, perplexity = 10, theta = 0.5, dims = 2)
plot(tmp.tsne$Y, main="Preeclampsia")
dev.off()
tmp.tsne <- Rtsne(t(tmp), check_duplicates = FALSE, pca = TRUE, perplexity = 10, theta = 0.5, dims = 2)
plot(tmp.tsne$Y, pch = 19, cex = 2, col = mypal[1:length(unique(pe.pd$Country))][pe.pd$Country])
text(tmp.tsne$Y, labels = pe.pd$GroupShort, col=c("black","darkred")[pe.pd$Group], cex=0.8)
plotMDS(tmp, top = nrow(tmp), labels = pe.pd$GroupShort, font=2,
        col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country])
## Look at samples relationship after ComBat
modcombat = model.matrix(~as.factor(Group), data=pe.pd)
pe.exprs.cb = ComBat(dat=pe.exprs, batch=factor(pe.pd$Dataset), 
                     mod=modcombat, par.prior=TRUE, prior.plots=TRUE)
dev.off(); rm(modcombat)
par(mfrow=c(2,1))
boxplot(pe.exprs, col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country], pch=".", main="Before")
boxplot(pe.exprs.cb, col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country], pch=".", main="After")
par(mfrow=c(1,2))
plotMDS(pe.exprs, top = nrow(pe.exprs), labels = pe.pd$GroupShort, font=2,
        col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country], main="Before")
plotMDS(pe.exprs.cb, top = nrow(pe.exprs.cb), labels = pe.pd$GroupShort, font=2,
        col=mypal[1:length(unique(pe.pd$Country))][pe.pd$Country], main="After")

## Differential expression analysis (DEA) of PE versus CTRL using limma
tmp.exprs <- pe.exprs
tmp.pd <- pe.pd
tmp.pd$Tissue <- factor(gsub("\\ ", "\\.", tmp.pd$Tissue))
by(tmp.pd$Group, tmp.pd$Dataset, table)
by(tmp.pd$Tissue, tmp.pd$Group, table)

## GLM design
design <- model.matrix(~0 + Group + Dataset, data = tmp.pd)
design
## Use array weights
aw <- arrayWeights(tmp.exprs, design)
fit <- lmFit(tmp.exprs, design, weights = aw)
contr.matrix <- makeContrasts(PreVsCtrl = GroupPreeclampsia - GroupControl, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contr.matrix), robust=TRUE)
contr.fit$genes <- pe.fd[rownames(contr.fit),]
summary(decideTests(contr.fit, p.value=0.1))
topTable(contr.fit, adjust = "BH")
write.table(topTable(contr.fit, adjust="BH", number = nrow(tmp.exprs)), 
            file = "DEA_Results_AW_PEvsCTRL.txt", row.names=FALSE, quote=FALSE, sep="\t")
tmp <- topTable(contr.fit, adjust = "BH", p.value = 0.1, number=nrow(tmp.exprs))
tmp <- tmp[tmp$adj.P.Val < 0.1, ]

## MDS plot using DEG signature
plotMDS(tmp.exprs[rownames(tmp),], labels = tmp.pd$GroupShort, 
        col=mypal[1:length(unique(pe.pd$Country))][factor(pe.pd$Country)])

## Plots for top differentially expressed genes 
myprobe <- rownames(tmp)[which.max(tmp$logFC)]
mygroup <- factor(tmp.pd$Group)
mybatch <- factor(tmp.pd$Dataset)
boxplot(tmp.exprs[myprobe,]~mygroup*mybatch, col = mypal2[c(1,4,2,6)][1:length(levels(mygroup))], main = pe.fd[myprobe, 'hgnc_symbol'])


#### Reanalysis of IUGR samples ####

iugr.pd$GroupShort <- iugr.pd$Group
levels(iugr.pd$GroupShort) <- c("C", "GR", "PE")
## Exclude Preeclampsia samples
iugr.exprs <- iugr.exprs[, iugr.pd$Group != "Preeclampsia"]
iugr.pd <- iugr.pd[iugr.pd$Group != "Preeclampsia", ]
iugr.pd$Group <- factor(iugr.pd$Group)
iugr.pd$GroupShort <- factor(iugr.pd$GroupShort)
## MDS plots
tmp <- iugr.exprs
pd <- iugr.pd
plotMDS(tmp,
        top = nrow(tmp), labels = pd$GroupShort,
        col = mypal[2:1][pd$Country])
legend("top", legend = levels(pd$Country), pch = 15, col = mypal[2:1])
plotMDS(tmp,
        top = nrow(tmp), labels = pd$Dataset,
        col = mypal[1:length(unique(pd$Country))][pd$Country])
plotMDS(tmp,
        top = nrow(tmp), labels = pd$GroupShort, pch = 19,
        col = mypal[1:length(unique(pd$Tissue))][pd$Tissue])
legend("bottomleft", legend = levels(pd$Tissue), pch = 15, col = mypal[1:length(unique(pd$Tissue))])
## Look at samples relationship after ComBat
modcombat <- model.matrix(~ as.factor(Group), data = iugr.pd)
iugr.exprs.cb <- ComBat(
  dat = iugr.exprs, batch = factor(iugr.pd$Dataset),
  mod = modcombat, par.prior = TRUE, prior.plots = TRUE)
dev.off()
rm(modcombat)
par(mfrow = c(2, 1))
boxplot(iugr.exprs, col = mypal[1:length(unique(pd$Country))][pd$Country], pch = ".", main = "Before")
boxplot(iugr.exprs.cb, col = mypal[1:length(unique(pd$Country))][pd$Country], pch = ".", main = "After")
par(mfrow = c(1, 2))
plotMDS(iugr.exprs,
        top = nrow(iugr.exprs), labels = pd$GroupShort, font = 2,
        col = mypal[1:length(unique(pd$Country))][pd$Country], main = "Before")
plotMDS(iugr.exprs.cb,
        top = nrow(iugr.exprs.cb), labels = pd$GroupShort, font = 2,
        col = mypal[1:length(unique(pd$Country))][pd$Country], main = "After")

## DEA of IUGR versus CTRL using limma
## Using blocking design with array weights
tmp <- iugr.exprs
pd <- iugr.pd
design <- model.matrix(~ 0 + Group + Dataset, data = pd)
design
aw <- arrayWeights(tmp, design)
names(aw) <- colnames(tmp)
fit <- lmFit(tmp, design, weights = aw)
contr.matrix <- makeContrasts(IugrVsCtrl = GroupIUGR - GroupControl, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contr.matrix))
contr.fit$genes <- iugr.fd[rownames(contr.fit), ]
summary(decideTests(contr.fit, p.value = 0.1))
topTable(contr.fit, adjust = "BH")
write.table(topTable(contr.fit, adjust = "BH", number = nrow(tmp)),
            file = "DEA_Results_AW_IUGRvsCTRL.txt", row.names = FALSE, quote = FALSE, sep = "\t")

## Plots for specific probes
myprobe <- "NM_000230"
mygroup <- factor(pd$Group)
mybatch <- factor(pd$Dataset)
boxplot(tmp[myprobe, ] ~ mygroup * mybatch, col = mypal2[c(1, 4, 2, 6)][1:length(levels(mygroup))])


## Reanalysis of PB samples

pb.pd$GroupShort <- pb.pd$Group
levels(pb.pd$Group)
levels(pb.pd$GroupShort) <- c("C", "PB")
by(pb.pd$Group, pb.pd$Dataset, table)

## MDS plots
tmp <- pb.exprs
pd <- pb.pd
plotMDS(tmp,
        top = nrow(tmp), labels = pd$GroupShort,
        col = mypal[1:2][pd$Country])
legend("bottomright", legend = levels(pd$Country), pch = 15, col = mypal[1:2])
plotMDS(tmp,
        top = nrow(tmp), labels = pd$Dataset,
        col = mypal[1:length(unique(pd$Country))][pd$Country])
plotMDS(tmp,
        top = nrow(tmp), labels = pd$GroupShort, pch = 19,
        col = mypal[1:length(unique(pd$Tissue))][pd$Tissue])
legend("bottomright", legend = levels(pd$Tissue), pch = 15, col = mypal[1:length(unique(pd$Tissue))])
plotMDS(tmp,
        top = nrow(tmp), labels = pd$GroupShort, pch = 19,
        col = mypal[1:length(unique(pd$Labor.Induction))][pd$Labor.Induction])
legend("bottomright", legend = levels(pd$Labor.Induction), pch = 15, col = mypal[1:length(unique(pd$Labor.Induction))])
## Look at samples relationship after ComBat
modcombat <- model.matrix(~ as.factor(Group), data = pb.pd)
pb.exprs.cb <- ComBat(
  dat = pb.exprs, batch = factor(as.character(svobj$sv > 0)),
  mod = modcombat, par.prior = TRUE, prior.plots = TRUE)
rm(modcombat)
par(mfrow = c(2, 1))
boxplot(pb.exprs, col = mypal[1:length(unique(pd$Country))][pd$Country], pch = ".", main = "Before")
boxplot(pb.exprs.cb, col = mypal[1:length(unique(pd$Country))][pd$Country], pch = ".", main = "After")
par(mfrow = c(1, 2))
plotMDS(pb.exprs,
        top = nrow(pb.exprs), labels = pd$GroupShort, font = 2,
        col = mypal[1:length(unique(pd$Country))][pd$Country], main = "Before")
plotMDS(pb.exprs.cb,
        top = nrow(pb.exprs.cb), labels = pd$GroupShort, font = 2,
        col = mypal[1:length(unique(pd$Country))][pd$Country], main = "After")

## DEA of PB versus CTRL using limma

## Pick all tissues
tmp <- pb.exprs
pd <- pb.pd
pd$Tissue <- factor(gsub("\\ ", "\\.", pd$Tissue))

## Pick only Placenta
tmp <- pb.exprs[, pb.pd$Tissue %in% "Placenta"]
pd <- pb.pd[pb.pd$Tissue %in% "Placenta", ]
pd$Tissue <- factor(pd$Tissue)

## Pick only Decidua
tmp <- pb.exprs[, pb.pd$Tissue %in% "Decidua basalis"]
pd <- pb.pd[pb.pd$Tissue %in% "Decidua basalis", ]
pd$Tissue <- factor(pd$Tissue)

## Design model
design <- model.matrix(~ 0 + Group, data = pd)
## Array weights
aw <- arrayWeights(tmp, design)
fit <- lmFit(tmp, design, weights = aw)
contr.matrix <- makeContrasts(PbVsCtrl = GroupPB - GroupControl, levels = design)
contr.fit <- eBayes(contrasts.fit(fit, contr.matrix))
contr.fit$genes <- pb.fd[rownames(contr.fit), ]
summary(decideTests(contr.fit, p.value = 0.1))
topTable(contr.fit, adjust = "BH")
write.table(topTable(contr.fit, adjust = "BH", number = nrow(tmp)),
            file = "DEA_Results_AW_PBvsCTRL.txt", row.names = FALSE, quote = FALSE, sep = "\t")

## Plots for specific probes
myprobe <- "ENSG00000110002"
mygroup <- factor(pd$Group)
mybatch <- factor(pd$Dataset)
boxplot(tmp[myprobe, ] ~ mygroup,
        col = mypal2[c(1, 4, 2, 6)][1:length(levels(mygroup))],
        main = paste0("Expression level of ", pb.fd[myprobe, "hgnc_symbol"], " gene"))


#### Processed data integration and exploratory analysis ####

all.probes <- Reduce(intersect, list(pe.fd$ensembl_gene_id, pb.fd$ensembl_gene_id, iugr.fd$ensembl_gene_id))
all.fd <- pe.fd[all.probes, ]
all.pd <- rbind(data.frame(pe.pd[, c(1:5, 10:11)], DS.Type = "PE"), 
                data.frame(pb.pd, DS.Type = "PB"), 
                data.frame(iugr.pd[, c(1:5, 8:9)], DS.Type = "IUGR"))
all.exprs <- cbind(pe.exprs[all.probes, ], pb.exprs[all.probes, ], iugr.exprs[all.probes, ])
# Exclude duplicated samples
all.pd <- all.pd[!(duplicated(all.pd$Sample)),]
all.exprs <- all.exprs[, all.pd$Sample]
tmp.pd <- all.pd
tmp.exprs <- all.exprs
## Average intensity level
plot(density(colMeans(tmp.exprs)))
## Intensity boxplots
my.group <- as.character(tmp.pd$Dataset)
boxplot(tmp.exprs, col = mypal3[factor(my.group)], pch=".", 
        main="Processed signal intensity")
legend("topright", legend = levels(factor(my.group)), pch = 15,
       col = mypal3[1:length(levels(factor(my.group)))], inset = c(-0.1, 0))
## MDS plot
png("exprs_mds.png", width = 3440, height = 1440, res = 150, units = "px",
    type = "cairo-png")
plotMDS(tmp.exprs, labels = as.character(tmp.pd$GroupShort), font = 2,
        col = mypal3[factor(my.group)], main = "Processed data MDS")
dev.off()

## Plot PCA with grouping using ggplot2
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggforce)
library(grDevices)
tmp.pca <- prcomp(t(tmp.exprs))
tmp.explained <- ((tmp.pca$sdev^2)/(sum(tmp.pca$sdev^2)))*100
tmp_plot <- tmp.pd %>% 
  mutate(Tissue = recode_factor(Tissue, "Chorionic villi" = "Chorion"),
         "Группа" = recode_factor(
           Group, "Control" = "Контроль", "Preeclampsia" = "Преэклампсия", 
           "PB" = "Преждевременные роды", 
           "IUGR" = "Задержка роста плода"), 
         "Ткань" = recode_factor(
           Tissue, "Chorion" = "Хорион", "Placenta" = "Плацента"), 
         PC1 = tmp.pca$x[,1], PC2 = tmp.pca$x[,2])
p <- tmp_plot %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_mark_ellipse(aes(fill = Dataset, label = Dataset), label.fontface = "plain",
                    show.legend = FALSE, alpha = 0.08, colour = "darkgrey", linetype = "dotted") +
  scale_shape_manual(values = c(0:3)) +
  geom_point(aes(color = Tissue, shape = Group), size = 2, stroke = 1) + 
  scale_color_manual(values = mypal3) +
  labs(title = NULL,
       x = paste0("PC1 (", round(tmp.explained[1], 2), "%)"),
       y = paste0("PC2 (", round(tmp.explained[2], 2), "%)")) +
  theme_bw() +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))
ggsave(plot = p, filename = "pca_exprs_datasets.pdf", 
       dpi = 300, width = 14, height = 9)
p2 <- tmp_plot %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_mark_ellipse(aes(fill = Dataset, label = Dataset), label.fontface = "plain",
                    show.legend = FALSE, alpha = 0.08, colour = "darkgrey", linetype = "dotted") +
  scale_shape_manual(values = c(0:3)) +
  geom_point(aes(color = Ткань, shape = Группа), size = 2, stroke = 1) + 
  scale_color_manual(values = mypal3) +
  labs(title = NULL,
       x = paste0("Компонента 1 (", round(tmp.explained[1], 0), "% дисперсии)"),
       y = paste0("Компонента 2 (", round(tmp.explained[2], 0), "% дисперсии)")) +
  theme_bw() +
  theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14, face = "bold"))
ggsave(plot = p2, filename = "pca_exprs_datasets_ru.pdf", 
       device = cairo_pdf, dpi = 300, width = 14, height = 9)
