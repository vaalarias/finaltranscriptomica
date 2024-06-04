library(Rsubread)
library(dplyr)
filess = c("RNA-hipocampo-5xFAD-A.bam", "RNA-hipocampo-5xFAD-B.bam", "RNA-hipocampo-WT-A.bam", "RNA-hipocampo-WT-B.bam")
ftchisat <- featureCounts(files=filess, annot.ext="mm10-ensembl_99-genes.gtf", isGTFAnnotationFile = T, useMetaFeatures = T, minMQS = 10, largestOverlap = T, isPairedEnd = T, requireBothEndsMapped = T, nthreads = 5)
glimpse(ftchisat)

### EDGER

install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)

counts <- ftchisat$counts

# Crear un vector con los nombres de las muestras
sample_names <- colnames(counts)

# Definir los grupos experimentales (ajusta esto según tu diseño experimental)
# Por ejemplo, "5xFAD" para las muestras de ratones 5xFAD y "WT" para las muestras de tipo salvaje
group <- factor(c("5xFAD", "5xFAD", "WT", "WT"))

genes=ftchisat$annotation

# Crear el objeto DGEList
dge <- DGEList(counts=counts, group=group, genes=genes[,c("GeneID","Length")])

# Visualizar el objeto DGEList para verificar
dge

# Filtrado de genes con baja expresión
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes=FALSE]

# Normalización de los datos
dge <- calcNormFactors(dge)

# Establecer el diseño experimental para el análisis diferencial
design <- model.matrix(~group)

# Estimación de la dispersión
dge <- estimateDisp(dge, design)

# Realización del análisis de expresión diferencial
fit <- glmFit(dge, design)
lrt <- glmLRT(fit)

# Ver los resultados del análisis
topTags(lrt)

# Guardar los resultados en un archivo
write.csv(topTags(lrt, n=Inf)$table, file="edgeR_results.csv")