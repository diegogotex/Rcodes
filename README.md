GOtools
================

### The Project

This project aims to help users to analyse and procuce good outputs for gene set enrichment analysis.

The enrichment process comes from other R packages like **ClusterProfiler** and **enrichR** and here we just work to give a better visualization of the results from these packages.

### Dataset

For an initial step lets select a data set to workthrough.

|               |  log2FoldChange|       padj|  exp|
|---------------|---------------:|----------:|----:|
| DDX11L1       |       -1.125993|  0.0495855|   -1|
| PLEKHN1       |        1.072670|  0.0358392|    1|
| HES4          |        1.868776|  0.0058709|    1|
| ISG15         |        4.598033|  0.0000000|    1|
| AGRN          |        2.391258|  0.0000000|    1|
| RP11-465B22.3 |        1.003380|  0.0104877|    1|
| RP1-140A9.1   |        1.545089|  0.0149719|    1|
| TP73          |        2.360830|  0.0000276|    1|
| AJAP1         |       -2.364369|  0.0000006|   -1|
| CHD5          |       -1.110567|  0.0356134|   -1|

### Enrichment

#### [**ClusterProfiler**](http://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)

Here i'm showing an exemple of an enhichment analysis using the function **enrichGO** from **ClusterProfiler**.

``` r
library(clusterProfiler)
library(org.Hs.eg.db)

GO_CP <- enrichGO(gene = DE_agudo$symbol,  #column from the dataframe with the symbols
                  OrgDb = org.Hs.eg.db,    #orgdb
                  keyType = "SYMBOL",      #type of data
                  ont = "BP")              #especify the GO ontology 
```

#### [**enrichR**](https://cran.r-project.org/web/packages/enrichR/index.html)

Here is an example of an enrichment analysis using the enrichR.

``` r
library(enrichR)

GO_ER <- enrichr(genes = DE_agudo$symbol, #column from the dataframe with the symbols
                 databases = "GO_Biological_Process_2018" #database to perform the enrichment analysis
                 )
```

    Uploading data to Enrichr... Done.
      Querying GO_Biological_Process_2018... Done.
    Parsing results... Done.

#### coolTable

This function recevies the output from enrichR or clusterProfiler and adds some information to the dataframe. it creates four columns called "up", "down", "up\_genes" and "down\_genes", related to \# of genes upregulates, \# of genes downregulated, a vector of genes up and downregulated, respectively. This function requires the specification of an enrichment type for the input object in **enrich.type**, choose **1 for clusterProfiler** input and **2 for enrichR** input.

``` r
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

                       ID
    GO:0006958 GO:0006958
    GO:0002455 GO:0002455
    GO:0006956 GO:0006956
    GO:0016064 GO:0016064
    GO:0019724 GO:0019724
    GO:0072376 GO:0072376
                                                                  Description
    GO:0006958                       complement activation, classical pathway
    GO:0002455 humoral immune response mediated by circulating immunoglobulin
    GO:0006956                                          complement activation
    GO:0016064                        immunoglobulin mediated immune response
    GO:0019724                                       B cell mediated immunity
    GO:0072376                                     protein activation cascade
               Term_size     p.adjust up down
    GO:0006958       137 7.451780e-69 86    0
    GO:0002455       149 1.813608e-67 87    1
    GO:0006956       170 2.923462e-58 86    0
    GO:0016064       212 4.152903e-55 89    3
    GO:0019724       213 5.547052e-55 89    3
    GO:0072376       194 1.283225e-53 87    0
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      up_genes
    GO:0006958                C1QA/C1QC/C1QB/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7
    GO:0002455           C1QA/C1QC/C1QB/EXO1/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7
    GO:0006956                C1QA/C1QC/C1QB/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7
    GO:0016064 C1QA/C1QC/C1QB/IL10/EXO1/IRF7/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7
    GO:0019724 C1QA/C1QC/C1QB/IL10/EXO1/IRF7/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7
    GO:0072376            C1QA/C1QC/C1QB/SERPING1/IGHA2/IGHG4/IGHG2/IGHA1/IGHG1/IGHG3/IGHM/IGHV6-1/IGHV1-3/IGHV2-5/IGHV3-7/IGHV3-11/IGHV3-13/IGHV3-15/IGHV1-18/IGHV3-20/IGHV3-21/IGHV3-23/IGHV1-24/IGHV2-26/IGHV4-28/IGHV3-30/IGHV3-33/IGHV4-34/IGHV4-39/IGHV3-43/IGHV3-48/IGHV3-49/IGHV5-51/IGHV3-53/IGHV1-58/IGHV4-59/IGHV4-61/IGHV3-64/IGHV3-66/IGHV1-69/IGHV2-70/IGHV1-69-2/IGHV3-72/IGHV3-73/IGHV3-74/IGKC/IGKV4-1/IGKV5-2/IGKV1-5/IGKV1-12/IGKV3-15/IGKV1-16/IGKV1-17/IGKV3-20/IGKV2-28/IGKV2-29/IGKV2-30/IGKV1-39/IGKV2-40/IGKV1D-39/IGKV1D-33/IGKV2D-30/IGKV2D-28/IGKV3D-20/IGKV1D-12/IGKV3D-11/IGLV6-57/IGLV1-51/IGLV1-47/IGLV1-44/IGLV7-43/IGLV1-40/IGLV3-27/IGLV3-25/IGLV2-23/IGLV3-21/IGLV3-19/IGLV2-14/IGLV2-11/IGLV2-8/IGLV3-1/IGLL5/IGLC2/IGLC3/IGLC6/IGLC7/GP9
                        down_genes
    GO:0006958                    
    GO:0002455            HLA-DQB1
    GO:0006956                    
    GO:0016064 TNFSF4/IL4/HLA-DQB1
    GO:0019724 TNFSF4/IL4/HLA-DQB1
    GO:0072376
