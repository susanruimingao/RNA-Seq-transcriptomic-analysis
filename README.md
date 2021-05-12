# RNA-Seq transcriptomic analysis

1) fastp: for raw reads QC and obtain trimmed reads if rquired:
```
for i in *_R1.fastq.gz ; do j=$(basename $i _R1.fastq.gz) ; fastp -i ${j}_R1.fastq.gz -I ${j}_R2.fastq.gz -o ${j}_trimmed_R1.fastq.gz -O ${j}_trimmed_R2.fastq.gz -h ${j}_fastp.htmp -j ${j}_fastp.json ; done
```
2) RSEM: reads mapping and expression calculation

I. Preparing Reference Sequences
```
rsem-prepare-reference Bm_atcc23344cDNA.fa Bm_atcc23344cDNA
```
II. Calculating Expression Values
```
for i in *_R1.fastq.gz ; do j=$(basename $i _R1.fastq.gz) ; echo "rsem-calculate-expression -p 12 --bowtie2 -paired-end ${j}_R1.fastq.gz ${j}_R2.fastq.gz Bm_atcc23344cDNA ${j}" ; done
```

III. data statistics learned by RSEM:
```
for i in *.genes.results ; do j=$(basename $i .genes.results) ; rsem-plot-model $j $j.pdf ; done
```

3) in R studio, install the required packages and load necessary packages
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("tximport")
BiocManager::install("tximportData")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("csaw")
BiocManager::install("mixOmics")

library(tximportData)
library(tximport)
library(limma)
library(csaw)
library(edgeR)
library(csaw)
library(RColorBrewer)
library(mixOmics)
library(magrittr)
library(tidyverse)

setwd("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/")

samples <- read.table(file.path("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/Bmallei_meta.txt"), header = TRUE)
samples

files <- file.path("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results", paste0(samples$isolate, ".genes.results"))
names(files) <- samples$isolate
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
cts <- txi.rsem$counts
head(cts)

write.table(cts, file = "counts.txt", sep = "\t")

counts <- txi.rsem$counts

d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
d0

cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

group <- as.factor(samples$description)

plotMDS(d, col = as.numeric(group))
dev.off()

mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)

tmp <- voom(d0, mm, plot = T)


fit <- lmFit(y, mm)
head(coef(fit))

contr <- makeContrasts(groupWTstd_23344 - groupDD3008, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)


top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)

length(which(top.table$adj.P.Val < 0.05))
filDEG <- length(which(top.table$adj.P.Val < 0.05))
view(filDEG)

top.table %>%
  filter(P.Value < 0.05) %>%
  write.table(file = "filDEG_WT23344_DD3008.txt", row.names = F, sep = "\t", quote = F)


View(top.table)

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "WT23344_DD3008.txt", row.names = F, sep = "\t", quote = F)


top.table %>%
  filter(P.Value < 0.05,
         abs(logFC) > 1) %>% dim()



rawCountTable <- read.table("counts.txt", header=TRUE, sep="\t", row.names=1)
sampleInfo <- read.table("Bmallei_meta.txt", header=TRUE, sep="\t", row.names=1)
head(rawCountTable)
head(sampleInfo)

degFull <- DGEList(rawCountTable, group=sampleInfo$description)
degFull

pseudoCounts <- log2(degFull$counts+1)
head(pseudoCounts)

hist(pseudoCounts[,"X23344.1"])
boxplot(pseudoCounts, col="gray", las=3)

par(mfrow=c(1,2))
avalues <- (pseudoCounts[,11] + pseudoCounts[,12])/2
mvalues <- (pseudoCounts[,11] - pseudoCounts[,12])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="WT")
abline(h=0, col="red")
avalues <- (pseudoCounts[,20] + pseudoCounts[,21])/2
mvalues <- (pseudoCounts[,20] - pseudoCounts[,21])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="mutant")
abline(h=0, col="red")

sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists

cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cim(sampleDists, color=cimColor, symkey=FALSE)

dev.off()
```

## Obtaining DEGs with complicated experimental design
```
FourCounts<- read.table("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/countsWT_3strains.txt",
                        header = T, sep="\t", row.names = 1)
View(FourCounts)
FourMeta <- read.table("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/targets.txt", header = T, sep="\t")
View(FourMeta)
Group<- as.factor(FourMeta$description)

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
row.names(design)<- row.names(FourMeta)

contrasts<- makeContrasts(
  WT_23344_vs_Bogor= WT_23344 - Bogor,
  WT_23344_vs_Muk= WT_23344 - Muk,
  WT_23344_vs_Zagreb= WT_23344 - Muk,
  WT_23344_vs_BogorMuk= WT_23344 - (Bogor+Muk)/2,
  WT_23344_vs_MukZagreb= WT_23344 - (Muk+Zagreb)/2,
  WT_23344_vs_BogorMukZagreb= WT_23344 - (Bogor+Muk+Zagreb)/3,
  etc...
  levels= design
)

FourCounts[FourCounts<0]
str(FourCounts)
#Fit model
d <- DGEList(counts = FourCounts, group = Group)
d <- calcNormFactors(d)
d<- estimateGLMCommonDisp(d, design)
d<- estimateGLMTagwiseDisp(d, design)
fit<- glmFit(d, design)


lrt<- glmLRT(fit, contrast= contrasts[, "WT_23344_vs_Bogor"])
detable<- topTags(lrt)
View(detable$table)

lrt<- glmLRT(fit, contrast= contrasts[, "WT_23344_vs_BogorMukZagreb"])
detable<- topTags(lrt)
View(lrt)

dea <- lrt$table
View(dea)
write.table(dea, file = "DEGs_WT_23344_vs_BogorMukZagreb.txt", row.names = T, sep = "\t", quote = F)
```

## Principal Component Analysis (PCA) using ggfortify and ggplot2
```
if(!require(devtools)) install.packages("devtools")
devtools::install_github("sinhrks/ggfortify")
install.packages("ggfortify")

library (dplyr)
library(tidyverse)
library(ggplot2)
library(ggfortify)

df <- read.table(file.path("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/counts.txt"))
samples <- read.table(file.path("/isilon/cfia-ottawa-fallowfield/users/gaoru/Bmallei_RNA-seq/genes_results/Bmallei_meta_old.txt"), header = TRUE)

df %>% 
  filter(rowSums(.) > 0) %>% 
t() %>% 
  as.data.frame() %>% 
  cbind.data.frame(samples) -> aa;
ncol(aa)
names(samples)
pca_res <- prcomp(aa[1:4995], scale. = TRUE)

autoplot(pca_res, data = aa, colour = 'description', label = TRUE)
```
