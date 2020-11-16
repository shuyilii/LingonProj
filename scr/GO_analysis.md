---
title: "LingonProj gene ontology analysis"
author:
   name: "Shuyi Li"
   email: shuyi.li@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "16 November, 2020"
output:
  html_document:
    keep_md: true
---



## Seperate genes based on different patterns

Input for run_select_pattern.sh: RSEM result (GeneMat_HFD_Lingon_LDF.Ebseqresults/GeneMat_HFD_Lingon_LDF.Ebseqresults_FDR_0.05.tab)
Output: Gene list with EnsemblID	and GeneSymbol, classified as diffent patterns
(Output folder: Matrix/pattern)


```bash
scr/run_select_pattern.sh
```

## GO analysis

### Load data

```r
pattern1 <- read.table("Matrix/pattern/pattern1.csv", header = TRUE)
pattern2 <- read.table("Matrix/pattern/pattern2.csv", header = TRUE)
pattern3 <- read.table("Matrix/pattern/pattern3.csv", header = TRUE)
pattern4 <- read.table("Matrix/pattern/pattern4.csv", header = TRUE)
pattern2.FDR0.05 <- read.table("Matrix/pattern/pattern2_FDR0.05.csv", header = TRUE)
pattern3.FDR0.05 <- read.table("Matrix/pattern/pattern3_FDR0.05.csv", header = TRUE)
pattern4.FDR0.05 <- read.table("Matrix/pattern/pattern4_FDR0.05.csv", header = TRUE)
```

### Go analysis using clusterProfiler
(clusterProfiler v3.16.1  For help: https://guangchuangyu.github.io/software/clusterProfiler)

Download genome wide annotation database for Mouse

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
```

Go analysis

```r
library(clusterProfiler)
library(org.Mm.eg.db)
ego<- function(pattern_data){
   ego_pattern <- enrichGO(gene = pattern_data[,1], OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = 'ALL',       
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

   str1 <- "GO_results/"
   str2 <- "_GO_enrich.csv"
   file_name <- paste(str1, deparse(substitute(pattern_data)), str2, sep = "")
   write.csv(as.data.frame(ego_pattern),file_name,row.names = F)
   return(ego_pattern)
}
```
GO analysis for gene assigned to pattern1

```r
ego_pattern1 <- ego(pattern1)
par(mfrow = c(1,2))
barplot(ego_pattern1,showCategory=20,drop=T)
dotplot(ego_pattern1,showCategory=20)
```
GO analysis for gene assigned to pattern2

```r
ego_pattern2 <- ego(pattern2)
par(mfrow = c(1,2))
barplot(ego_pattern2,showCategory=20,drop=T)
dotplot(ego_pattern2,showCategory=20)
```
GO analysis for gene assigned to pattern2.FDR0.05

```r
ego_pattern2.FDR0.05 <- ego(pattern2.FDR0.05)
par(mfrow = c(1,2))
barplot(ego_pattern2.FDR0.05,showCategory=20,drop=T)
dotplot(ego_pattern2.FDR0.05,showCategory=20)
```
GO analysis for gene assigned to pattern3

```r
ego_pattern3 <- ego(pattern3)
par(mfrow = c(1,2))
barplot(ego_pattern3,showCategory=20,drop=T)
dotplot(ego_pattern3,showCategory=20)
```
GO analysis for gene assigned to pattern3.FDR0.05

```r
ego_pattern3.FDR0.05 <- ego(pattern3.FDR0.05)
par(mfrow = c(1,2))
barplot(ego_pattern3.FDR0.05,showCategory=20,drop=T)
dotplot(ego_pattern3.FDR0.05,showCategory=20)
```

GO analysis for gene assigned to pattern4

```r
ego_pattern4 <- ego(pattern4)
par(mfrow = c(1,2))
barplot(ego_pattern4,showCategory=20,drop=T)
dotplot(ego_pattern4,showCategory=20)
```
GO analysis for gene assigned to pattern4.FDR0.05

```r
ego_pattern4.FDR0.05 <- ego(pattern4.FDR0.05)
par(mfrow = c(1,2))
barplot(ego_pattern4.FDR0.05,showCategory=20,drop=T)
dotplot(ego_pattern4.FDR0.05,showCategory=20)
```

