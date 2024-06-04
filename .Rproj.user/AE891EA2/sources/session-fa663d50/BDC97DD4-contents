### FEATURE COUNTS
## needed libraries
library(Rsubread)
filess = c("RNA-hipocampo-5xFAD-A.bam", "RNA-hipocampo-5xFAD-B.bam", "RNA-hipocampo-WT-A.bam", "RNA-hipocampo-WT-B.bam")
ftcbowtie <- featureCounts(files=filess, annot.ext="mm10-ensembl_99-genes.gtf", isGTFAnnotationFile = T, useMetaFeatures = T, minMQS = 10, largestOverlap = T, isPairedEnd = T, requireBothEndsMapped = T, nthreads = 5)
ftcbowtie
## File save
count_data <- ftcbowtie$counts
genes_data <- ftcbowtie$annotation
head(genes_data)
write.csv(count_data, file = "counts_data.csv")
write.csv(genes_data, file = "genes_data.csv")
### DEA EdgeR
## needed libraries
install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
library(pheatmap)
library(ggplot2)

## table counts
counts <- ftcbowtie$counts

## dataprep
sample_names <- colnames(counts)
group <- factor(c("FAD", "FAD", "WT", "WT"))
genes <- ftcbowtie$annotation

## DGEList
dge <- DGEList(counts=counts, group=group, genes=genes[,c("GeneID", "Length")])

## Low expression genes
keep <- rowSums(cpm(dge) > 1) >= 3
sum(keep)
dge <- dge[keep,]

## Data normalization
dge <- calcNormFactors(dge)


## PCAs
color_all = c("blue","blue","red","red")
color_samples=c("blue", "red")
plotMDS(dge, cex=0.8, col=color_all)
#dev.off()

## Principal components
PCA_log2CPM = prcomp(t(cpm(dge, log=T)), center=T, scale=T)
#summary( PCA_log2CPM )
PCA_perc_var = round ( ((PCA_log2CPM$sdev^2 / sum(PCA_log2CPM$sdev^2))*100), 1)
barplot(PCA_perc_var, names=colnames(PCA_log2CPM$x), xlab="Principal components", ylab="Percent Variation")



## Model matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
## Dispersion
dge <- estimateDisp(dge, design)
## Ajuste del modelo
fit <- glmFit(dge, design)

## Contrasts
contrast <- makeContrasts(FAD - WT, levels=design)

## DEA usando el contraste
lrt <- glmLRT(fit, contrast=contrast)

## Results
results <- topTags(lrt, n=Inf)
diff_genes <- results$table

## Guardar los resultados en un archivo CSV
write.csv(diff_genes, file="edgeR_results.csv")

## Ver los primeros genes diferencialmente expresados
head(diff_genes)

### Volcano Plot
diff_genes$Significant <- diff_genes$FDR < 0.05 & abs(diff_genes$logFC) > 1

## Plotting
ggplot(diff_genes, aes(x=logFC, y=-log10(FDR), color=Significant)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values=c("black", "red")) +
  labs(x="Log2 Fold Change", y="-Log10 FDR", title="Volcano Plot")

### HeatMap
# Genes diferencialmente expresados
significant_genes <- diff_genes[diff_genes$FDR < 0.05 & abs(diff_genes$logFC) > 1, ]
significant_gene_names <- rownames(significant_genes)

# Log2 CPM de genes diferencialmente expresados
log2_cpm_diff_genes <- cpm(dge, log=TRUE)[significant_gene_names, ]

# Heatmap
pheatmap(log2_cpm_diff_genes, scale="row", show_rownames=FALSE)


