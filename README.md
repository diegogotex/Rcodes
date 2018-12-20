---
title: "dieGOtools"
fontsize: 12pt
linkcolor: blue
mainfont: Arial
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(xtable)
```
### The Project

This project aims to help users to analyse and procuce good outputs for gene set enrichment analysis. 

The enrichment process comes from other R packages like **ClusterProfiler** and **enrichR** and here we just work to give a better visualization of the results from these packages.

### Dataset

For an initial step lets select a data set to workthrough.

```{r echo=F, results='asis', comment=NA, message=FALSE}
load("~/Desktop/teste.RData")
knitr::kable(DE_agudo[1:10, c(3,7,9)], caption = "Table 1. exemple of a table with DEG.")
```

### Enrichment




#### [**ClusterProfiler**](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

Here i'm showing an exemple of an enhichment analysis using the function **enrichGO** from **ClusterProfiler**.


```{r, comment=NA, message=FALSE}
library(clusterProfiler)
library(org.Hs.eg.db)

GO_CP <- enrichGO(gene = DE_agudo$symbol,  #column from the dataframe with the symbols
                  OrgDb = org.Hs.eg.db,    #orgdb
                  keyType = "SYMBOL",      #type of data
                  ont = "BP")              #especify the GO ontology 


```




#### [**enrichR**](https://cran.r-project.org/web/packages/enrichR/index.html)

Here is an example of an enrichment analysis using the enrichR.

```{r, comment=NA, message=FALSE}
library(enrichR)

GO_ER <- enrichr(genes = DE_agudo$symbol, #column from the dataframe with the symbols
                 databases = "GO_Biological_Process_2018" #database to perform the enrichment analysis
                 )


```





#### coolTable

This function recevies the output from enrichR or clusterProfiler and adds some information to the dataframe. it creates four columns called "up", "down", "up_genes" and "down_genes", related to # of genes upregulates, # of genes downregulated, a vector of genes up and downregulated, respectively.  This function requires the specification of an enrichment type for the input object in **enrich.type**, choose **1 for clusterProfiler** input and **2 for enrichR** input.

```{r, comment=NA, message=FALSE}
source("~/Dropbox/Rcodes/cooltable.R")


#creating a vector with upregulated SYMBOL
up_agudo <- subset(DE_agudo$symbol, DE_agudo$log2FoldChange > 0)
#creating a vector with downregulated SYMBOL
down_agudo <- subset(DE_agudo$symbol, DE_agudo$log2FoldChange < 0)


CT_GO_CP <- cool_table(enrichment.obj = GO_CP, #enrichment object from clusterProfiler
                       enrich.type = 1,        #selecting clusterProfiler object as input
                       up = up_agudo,          #vector with upregulated genes
                       down = down_agudo       #vector with downregulated genes
                       )

head(CT_GO_CP, 6)
```

