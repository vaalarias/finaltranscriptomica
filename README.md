### Valentina Janet Arias Ojeda                                                     *Transcriptómica*

En este proyecto el objetivo es realizar un análisis transcriptómico de acuerdo al pipeline que se ha estudiado a lo largo del curso. Para este fin, se utilizó un dataset de la base de datos recount3, proveniente del estudio **“The BCR-ABL1 Inhibitors Imatinib and Ponatinib Decrease Plasma Cholesterol and Atherosclerosis, and Nilotinib and Ponatinib Activate Coagulation in a Translational Mouse Model”**

Los inhibidores de la tirosina quinasa BCR-ABL1, como Imatinib, Ponatinib y Nilotinib, son terapias fundamentales para el tratamiento de la leucemia mieloide crónica (LMC) y otras neoplasias hematológicas. Estas terapias dirigidas funcionan al bloquear la actividad de la proteína BCR-ABL1, una enzima que promueve la proliferación descontrolada de células cancerosas. Aunque estos fármacos han transformado el tratamiento de la LMC, su uso a largo plazo ha revelado una variedad de efectos secundarios que afectan significativamente la calidad de vida de los pacientes .

La proteína BCR-ABL1 es el resultado de una translocación cromosómica que crea una fusión entre los genes BCR y ABL1. Esta proteína fusionada tiene actividad tirosina quinasa constitutivamente activa, lo que promueve la proliferación celular y la inhibición de la apoptosis, características principales de la LMC . Los inhibidores de BCR-ABL1, como Imatinib, fueron los primeros en introducirse en la práctica clínica y han demostrado ser altamente efectivos en la inducción de remisiones hematológicas y citogenéticas completas.

El análisis detallado de los efectos secundarios de estos inhibidores es crucial para mejorar la terapia dirigida y minimizar los riesgos para los pacientes. Entender cómo estos fármacos alteran la expresión génica en tejidos relevantes puede proporcionar información sobre los mecanismos moleculares subyacentes a estos efectos secundarios. Este estudio utiliza datos de RNA-seq para investigar los cambios en la expresión génica inducidos por los inhibidores de BCR-ABL1 en un modelo de ratón. La integración de estos datos puede ayudar a identificar nuevas estrategias para manejar los efectos adversos y mejorar la seguridad del tratamiento a largo plazo .

### Librerías

Para la manipulación de los datos, hacer el análisis estadístico y crear las visualizaciones se utilizaron las siguientes librerías

```r
library("recount3")
library("edgeR")
library("limma")
library("ggplot2")
library("pheatmap")
library("dplyr")
```

## Metodología

**Preparación de los datos**

Se descargaron y procesaron los datos de RNA-seq utilizando la biblioteca `recount3` en R. El proceso de preparación de datos incluye varios pasos cruciales:

**Cálculo de Conteos de Lecturas**:
El primer paso es calcular los conteos de lecturas, que representan el número de veces que cada fragmento de RNA es leído durante el experimento de secuenciación. Esto se hace para cuantificar la expresión génica en cada muestra.

```r
## GetCounts
assay(rse_gene_SRP117739, "counts") <- compute_read_counts(rse_gene_SRP117739)
rse_gene_SRP117739$sra.sample_attributes
```

**Expansión de Atributos SRA**:
Los datos descargados del SRA contienen atributos experimentales codificados que necesitan ser descompuestos y expandidos para un análisis más detallado. Estos atributos incluyen información como el tipo de tratamiento, la dosis y el tejido de origen.

```r
## expand data and incorporate them to the dataframe
rse_gene_SRP117739 <- expand_sra_attributes(rse_gene_SRP117739)
colData(rse_gene_SRP117739)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP117739)))
]
```

**Conversión de Atributos a Factores**:
Para el análisis estadístico, los atributos categóricos se convierten en factores. Esto permite el uso adecuado de modelos lineales y otras técnicas estadísticas que requieren variables categóricas.

**Filtrado de Genes y Muestras**

Para asegurar la calidad y relevancia de los datos, se aplicaron varios filtros:

**Proporción de Genes Asignados**:
Se calculó la proporción de lecturas que se asignaron a genes específicos, en comparación con el total de lecturas. Las muestras con una baja proporción de lecturas asignadas pueden indicar problemas técnicos y se excluyen del análisis para evitar sesgos.

```r
rse_gene_SRP117739$assigned_gene_prop <- rse_gene_SRP117739$recount_qc.gene_fc_count_all.assigned / rse_gene_SRP117739$recount_qc.gene_fc_count_all.total
```

**Visualización de la Proporción de Genes Asignados**:
Se utilizó un histograma para visualizar la distribución de la proporción de genes asignados. Se estableció un umbral para excluir muestras de baja calidad.

```r
rse_gene_SRP117739 <- rse_gene_SRP117739[, rse_gene_SRP117739$assigned_gene_prop > 0.71]
```

**Filtrado Basado en la Expresión Génica Media**:
Para cada gen, se calculó la media de las lecturas a través de todas las muestras. Los genes con una baja expresión media se eliminaron para reducir el ruido y centrarse en genes biológicamente relevantes.

```r
gene_means <- rowMeans(assay(rse_gene_SRP117739, "counts"))
```

**Normalización y Análisis Diferencial**

La normalización y el análisis diferencial de expresión son cruciales para identificar genes que cambian significativamente entre las condiciones experimentales:

**Normalización**:
La normalización ajusta las diferencias en la profundidad de secuenciación entre muestras. Esto se hace usando el método TMM (Trimmed Mean of M-values), que es robusto frente a genes altamente expresados y outliers.

```r
## DGEList
dge <- DGEList(
  counts = assay(rse_gene_SRP117739, "counts"),
  genes = rowData(rse_gene_SRP117739)
)
dge
## Normalization
dge <- calcNormFactors(dge)
```

**Análisis de Expresión Diferencial**:

- **Modelo de Diseño**:
Se construye un modelo de diseño que incluye términos para las diferentes condiciones experimentales (tratamiento con distintos inhibidores). Este modelo permite estimar el efecto de cada tratamiento en la expresión génica.

```r
design <- model.matrix(~0+ sra_attribute.treatment_compound + assigned_gene_prop, data = colData(rse_gene_SRP117739))
```

- **Transformación Voom**:
La función `voom` transforma los datos de conteo en logaritmos base 2 de conteos por millón (CPM), ajustando la varianza para que sea aproximadamente constante en todo el rango de medias de expresión. Esto es necesario para aplicar modelos lineales y obtener estadísticas fiables.
    
    ![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/803b6be9-b008-4c17-ae0e-78403619a2c5/c382531c-ab4d-4a02-8da0-86027c89cffd/Untitled.png)
    
- **Modelo Lineal y Pruebas Estadísticas**:
Se ajusta un modelo lineal a los datos transformados usando `lmFit`. Luego, se aplica el método `eBayes` para obtener estadísticas bayesianas que mejoran la precisión de las estimaciones y reducen la tasa de falsos positivos.
    
    ![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/803b6be9-b008-4c17-ae0e-78403619a2c5/514df12e-4889-40f9-8366-60f88af90822/Untitled.png)
    
- **Obtención de Resultados**:
Se generan tablas con los genes más significativamente diferentes entre las condiciones, ordenados por p-valor ajustado. Esto permite identificar los genes que muestran cambios más robustos y significativos en respuesta a los tratamientos.

## Resultados

![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/803b6be9-b008-4c17-ae0e-78403619a2c5/c2409a65-7433-4fdc-8d8a-023bdaba0297/Untitled.png)

![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/803b6be9-b008-4c17-ae0e-78403619a2c5/03a04dbe-86c1-4e98-8e0a-3d7cc9f40e2c/Untitled.png)

![Untitled](https://prod-files-secure.s3.us-west-2.amazonaws.com/803b6be9-b008-4c17-ae0e-78403619a2c5/5cd3230c-6dc8-4ce2-b302-eb39aa8f12ce/Untitled.png)

El análisis del heatmap indica variaciones significativas en la expresión génica entre los tratamientos con Imatinib, Ponatinib y Nilotinib. Se observa una regulación diferencial de genes asociados con funciones biológicas clave, como la adhesión celular (Integrina β2 y α4, Icam1, Vcam1), la activación de macrófagos (Cd40, Tnfrsf14, Nfκb), la regulación lipídica (Apoa1, Apoa2, Apoc2, Apoc4), los procesos inflamatorios (Cd40, Nfκb, Tnfrsf14, Icam1, Vcam1) y la modulación de la matriz extracelular (Col1a2, Col3a1, Col1a1, Tgf-β). Estos genes son conocidos por desempeñar un papel fundamental en la señalización de aterosclerosis y otros procesos cardiovasculares. 

La mayor regulación positiva de estos genes por parte del Imatinib sugiere un impacto funcional más pronunciado en la mitigación de los riesgos cardiovasculares asociados con el tratamiento. Por otro lado, los efectos observados con Ponatinib y Nilotinib, aunque presentes, son menos marcados en comparación con Imatinib. Esta asociación entre la expresión génica y las funciones biológicas nos permite deducir el posible impacto funcional de cada tratamiento y resaltar la importancia de considerar estos aspectos en la selección y administración de terapias dirigidas para pacientes con leucemia mieloide crónica.

### Conclusiones

Este análisis reforzó los hallazgos reportados en el artículo publicado mediante estableciendo una base estadística más sólida a los resultados de cómo los inhibidores de BCR-ABL1 afectan la expresión génica relacionada con el metabolismo del colesterol y la coagulación. Los resultados sugieren que mientras algunos inhibidores pueden tener efectos beneficiosos sobre el colesterol y la aterosclerosis, otros pueden aumentar el riesgo de coagulación. Estos hallazgos son importantes para guiar futuras investigaciones y el desarrollo de estrategias terapéuticas para minimizar los efectos secundarios en pacientes tratados con estos fármacos.

**Referencias y reproducibilidad**

- repositorio en github: [vaalarias/finaltranscriptomica (github.com)](https://github.com/vaalarias/finaltranscriptomica)
- Pouwer, M. G., Pieterman, E. J., Verschuren, L., Caspers, M. P. M., Kluft, C., Garcia, R. A., Aman, J., Jukema, J. W., & Princen, H. M. G. (2018). The BCR-ABL1 inhibitors imatinib and ponatinib decrease plasma cholesterol and atherosclerosis, and nilotinib and ponatinib activate coagulation in a translational mouse model. *Frontiers in Cardiovascular Medicine*, *5*. https://doi.org/10.3389/fcvm.2018.00055
- Hochhaus, A., Baccarani, M., Silver, R. T., Schiffer, C., Apperley, J. F., Cervantes, F., ... & Simonsson, B. (2016). European LeukemiaNet 2020 recommendations for treating chronic myeloid leukemia. Leukemia, 34(4), 966-984.
- Druker, B. J., Tamura, S., Buchdunger, E., Ohno, S., Segal, G. M., Fanning, S., ... & Lydon, N. B. (1996). Effects of a selective inhibitor of the Abl tyrosine kinase on the growth of Bcr-Abl positive cells. Nature medicine, 2(5), 561-566.
- Robinson MD, Oshlack A: A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol.* 2010; **11**(3): R25
