### DEA con DESeq2
## Libraries needed
BiocManager::install("DESeq2")
BiocManager::install("airway")
BiocManager::install("tximport")
BiocManager::install("tidyverse")
library(DESeq2)
library(tidyverse)
library(airway)
library(tximport)
library(ggplot2)
library(pheatmap)

## DataPrep
count_data <- read.csv("counts_data.csv", row.names = 1)
#head(count_data)
sample_names <- colnames(counts)
col_data <- data.frame(
  sample_name = colnames(count_data),
  group = c("FAD", "FAD", "WT", "WT")
)
## DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = col_data,
                              design = ~ group)
#dds
## Filtering low counts
keep<-rowSums(counts(dds)) >=10
dds<- dds[keep,]
#dds
## Differential Exp Analysis
dds <- DESeq(dds)
res <- results(dds)
## Adjust pvalue
res0.01<- results(dds,alpha=0.01)

## Save results
write.csv(as.data.frame(res0.01), file = "deseq2_results.csv")

## PCA
rld <- rlog(dds, blind = FALSE)

## PCA Plot
pcaData <- plotPCA(rld, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()

## Barplot for PC
pca <- prcomp(t(assay(rld)))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
barplot(100 * percentVar, xlab = "Principal Component", ylab = "Percentage of Variance Explained")

## Volcano Plot
res$Significant <- res$padj < 0.05 & abs(res$log2FoldChange) > 1
ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 adjusted p-value", title = "Volcano Plot")

## Heatmap 
sig_genes <- rownames(res)[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)]
log2_cpm <- assay(rld)[sig_genes, ]
pheatmap(log2_cpm, scale = "row", show_rownames = FALSE)
