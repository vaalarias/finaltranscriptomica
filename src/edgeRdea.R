## needed libraries
library("recount3")
library("edgeR")
library("limma")
library("ggplot2")
library("pheatmap")
library("dplyr")
## Study: SRP117739
### DATA LOAD
recount3::create_rse_manual(
  project = "SRP117739",
  project_home = "data_sources/sra",
  organism = "mouse",
  annotation = "gencode_v23",
  type = "gene"
)
## Available projects
mouse_projects <- available_projects(organism = "mouse")
## Subset
rse_gene_SRP117739 <- create_rse(
  subset(
    mouse_projects,
    project == "SRP117739" & project_type == "data_sources"
  )
)

## GetCounts
assay(rse_gene_SRP117739, "counts") <- compute_read_counts(rse_gene_SRP117739)
rse_gene_SRP117739$sra.sample_attributes

### DATA PROCESSING

## expand data and incorporate them to the dataframe
rse_gene_SRP117739 <- expand_sra_attributes(rse_gene_SRP117739)

colData(rse_gene_SRP117739)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP117739)))
]
rse_gene_SRP117739$sra.sample_attributes
## Datatype to factor due to characters
rse_gene_SRP117739$sra_attribute.source_name <- factor(rse_gene_SRP117739$sra_attribute.source_name)
rse_gene_SRP117739$sra_attribute.tissue <- factor(rse_gene_SRP117739$sra_attribute.tissue)
rse_gene_SRP117739$sra_attribute.sra_attribute.treatment_compound <- factor(rse_gene_SRP117739$sra_attribute.treatment_compound)
rse_gene_SRP117739$sra_attribute.treatment_dose <- factor(rse_gene_SRP117739$sra_attribute.treatment_dose)
rse_gene_SRP117739$sra_attribute.treatment_short_descript <- factor(rse_gene_SRP117739$sra_attribute.treatment_short_descript)
rse_gene_SRP117739$sra_attribute.treatment_unit <- factor(rse_gene_SRP117739$sra_attribute.treatment_unit)
### GENE FILTERING

## reads proportion
rse_gene_SRP117739$assigned_gene_prop <- rse_gene_SRP117739$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP117739$recount_qc.gene_fc_count_all.total

## Summary by tissue group
with(colData(rse_gene_SRP117739), tapply(assigned_gene_prop, sra_attribute.treatment_compound, summary))

## Dots
with(colData(rse_gene_SRP117739), plot(assigned_gene_prop, sra_attribute.treatment_short_descript))
abline(v = 0.71, col = "deepskyblue4", lwd = 4, lty = "solid")
## Histrogram
hist(rse_gene_SRP117739$assigned_gene_prop, col="gray")
abline(v=0.70,col="purple", lwd=7, lty = "dashed")

## SAFE
rse_gene_SRP117739_copia<-rse_gene_SRP117739

## Get rid of samples <  0.70
table(rse_gene_SRP117739$assigned_gene_prop < 0.71)

## corrected object by read proportion
rse_gene_SRP117739 <- rse_gene_SRP117739[, rse_gene_SRP117739$assigned_gene_prop > 0.71]

## low significance genes
gene_means <- rowMeans(assay(rse_gene_SRP117739, "counts"))
summary(gene_means)

rse_gene_SRP117739 <- rse_gene_SRP117739[gene_means > 0.1, ]
#keep <- filterByExpr(rse_gene_SRP127581, group=group)
## Final Dimentions
dim(rse_gene_SRP117739)
## Retained gene percentage
round(nrow(rse_gene_SRP117739) / nrow(rse_gene_SRP117739_copia) * 100, 2)

####rse_gene_SRP127581<-rse_gene_SRP127581_copia

## DGEList
dge <- DGEList(
  counts = assay(rse_gene_SRP117739, "counts"),
  genes = rowData(rse_gene_SRP117739)
)
dge
## Normalization
dge <- calcNormFactors(dge)


### DEA
## boxplot 
ggplot(as.data.frame(colData(rse_gene_SRP117739)), aes(y = assigned_gene_prop, x = sra_attribute.treatment_compound )) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Treatment")
## Modelo Estadistico las lecturas de las cuentas a logaritmo base 2 de las cuentas por millón.
design <- model.matrix(~0+ sra_attribute.treatment_compound + assigned_gene_prop, data = colData(rse_gene_SRP117739))
design
vGene <- voom(dge, design, plot = TRUE)
# Aplicar eBayes para obtener estadísticas bayesianas
fit_bayes <- eBayes(lmFit(vGene))
plotMA(fit_bayes, coef = "sra_attribute.treatment_compoundCtrl", col = "steelblue")
results <- topTable(fit_bayes, coef = "sra_attribute.treatment_compoundCtrl", number = nrow(rse_gene_SRP117739), sort.by = "none")
results

# visualizacion
volcanoplot(fit_bayes, coef = "sra_attribute.treatment_compoundCtrl", highlight = 5, names = results$gene_name, col = "steelblue", hl.col="hotpink4")

exprs_heatmap <- vGene$E[rank(results$adj.P.Val) <= 30, ]
nombres <- rownames(results)
nombres
rownames(exprs_heatmap) <- results$gene_name[match(rownames(exprs_heatmap), nombres)]
df <- as.data.frame(colData(rse_gene_SRP117739)[, c("sra_attribute.treatment_compound","assigned_gene_prop" )])
colnames(df) <- c("treatment", "gene_prop")
df
pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_col = df
)
