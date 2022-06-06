-   [1. Load data into R as a Seurat
    object](#load-data-into-r-as-a-seurat-object)
-   [2. Verify data have been QC’d](#verify-data-have-been-qcd)
-   [3. Transform gene counts](#transform-gene-counts)

``` r
knitr::opts_chunk$set(echo = T, message = FALSE, warning = FALSE)
stdt<-date()

# Libraries
library(tidyverse)
```

    ## Warning in system("timedatectl", intern = TRUE): running command 'timedatectl'
    ## had status 1

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(Seurat)
```

    ## Attaching SeuratObject

    ## Attaching sp

``` r
library(cowplot)

# Paths
wd <- "/hpc/group/chsi-mic-2022"
datadirs <- list.files(path = file.path(wd, "data","Christian2021CellReports"), pattern = "_10x", full.names = T)
datadirs <- datadirs[!grepl(".h5", datadirs)]
intermeddir <- file.path(wd, "intermed")
outdir <- file.path("/work",Sys.info()[["user"]],"output")
outdir
```

    ## [1] "/work/mrl17/output"

**Goal of this workshop:** Learn how to load cellranger output into R
and prepare QC’d data for downstream analyses

**What’s covered in this workshop:** - Load data into R as a Seurat
object - Verify that the data have been QC’d - Transform gene counts

**Data source:** Christian et al. 2021 Cell Reports
(<https://doi.org/10.1016/j.celrep.2021.109118>)

Note that the Seurat - Guided Clustering Tutorial
(<https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>) includes
a lot of overlapping information.

**METHOD DETAILS, Single-cell RNA-seq analysis, page e4**

-   Pre-processing: “Raw short reads were demultiplexed, filtering and
    mapped to mouse genome GRCm38/mm10 using cellranger v2.02.”

-   QC: “The gene count matrices from cellranger were subjected to
    quality control, pre-processing and clustering using the R Seurat
    2.3.4 package (Butler et al., 2018). Low-quality cells that had less
    than 200 expressed genes and more than 5% mitochondrial genes were
    filtered out.”

-   Analysis: “Gene counts were scaled to total gene expression and
    percentage of mitochondrial genes with a scaling factor of 10,000,
    and then log-transformed. The high dimensional data for each sample
    were reduced by PCA and t-Distributed Stochastics Neighbor Embedding
    (tSNE). We used the FindCluster function to group clusters in each
    sample with a resolution of 0.6. Differential expressed genes (DEGs)
    were identified using the Wilcoxon rank-sum test.”

## 1. Load data into R as a Seurat object

We can use code like this to load the cell ranger output and store it in
a Seurat object.

``` r
readin_10x_to_SO <- function(sample, path, min.cells = 0, min.genes = 0){
  Seurat::Read10X(data.dir = path) %>%
    Seurat::CreateSeuratObject(project = sample, 
                               min.cells=min.cells, min.genes=min.genes) -> seu.obj
  return(seu.obj)
}
```

Load the four 10x datasets that correspond to T cells sampled in and
outside of the tumor. See Fig 2A for tSNE plots of each sample type.

-   Tem = effector memory in tumor
-   Tcm = central memory in tumor
-   Trm = resident memory in tumor
-   disTrm = resident memory T-cells outside of the tumor

``` r
sampleNames <- gsub("_10x","", basename(datadirs))
sampleNames
```

    ## [1] "disTrm" "Tcm"    "Tem"    "Trm"

``` r
seulist <- list()
for(i in 1:length(sampleNames)){
  seulist[[i]] <- readin_10x_to_SO(sample = sampleNames[i], path = datadirs[i])
}
names(seulist) <- sampleNames
seulist
```

    ## $disTrm
    ## An object of class Seurat 
    ## 12064 features across 863 samples within 1 assay 
    ## Active assay: RNA (12064 features, 0 variable features)
    ## 
    ## $Tcm
    ## An object of class Seurat 
    ## 12840 features across 4009 samples within 1 assay 
    ## Active assay: RNA (12840 features, 0 variable features)
    ## 
    ## $Tem
    ## An object of class Seurat 
    ## 12571 features across 1457 samples within 1 assay 
    ## Active assay: RNA (12571 features, 0 variable features)
    ## 
    ## $Trm
    ## An object of class Seurat 
    ## 11778 features across 459 samples within 1 assay 
    ## Active assay: RNA (11778 features, 0 variable features)

## 2. Verify data have been QC’d

*METHOD DETAILS, Single-cell RNA-seq analysis, page e4*: “Low-quality
cells that had less than 200 expressed genes and more than 5%
mitochondrial genes were filtered out.”

The rownames are gene symbols. Investigating the names a bit more, the
mitochondrial genes can be identified by the prefix “mt-”

``` r
tem <- seulist$Tem
curr.count <- tem@assays$RNA@counts
rownames(curr.count)[1:10]
```

    ##  [1] "Mrpl15"        "Lypla1"        "Tcea1"         "Atp6v1h"      
    ##  [5] "Rb1cc1"        "4732440D04Rik" "Pcmtd1"        "Gm26901"      
    ##  [9] "Rrs1"          "Mybl1"

``` r
rownames(curr.count)[grepl("^mt-", rownames(curr.count))] # these look like mitochondrial genes!
```

    ##  [1] "mt-Nd1"  "mt-Nd2"  "mt-Co1"  "mt-Co2"  "mt-Atp8" "mt-Atp6" "mt-Co3" 
    ##  [8] "mt-Nd3"  "mt-Nd4l" "mt-Nd4"  "mt-Nd5"  "mt-Nd6"  "mt-Cytb"

Add the mitochondrial percentage for all seurat objects in our list

``` r
add_mt <- function(so){
  so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
  return(so)
}

seulist %>%
  map(~add_mt(.x)) -> seulist.mt

seulist.mt %>%
  map(~head(.x@meta.data))
```

    ## $disTrm
    ##                  orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGAGTCCTCCT     disTrm       2459          963   3.538024
    ## AAACCTGCACTTAAGC     disTrm       3035         1038   3.228995
    ## AAACCTGTCTGATACG     disTrm       6070         1745   2.075783
    ## AAACGGGAGCGATCCC     disTrm       8625         1700   3.524638
    ## AAAGTAGCACATGTGT     disTrm       3969         1403   2.670698
    ## AAATGCCGTCAATACC     disTrm       2712          968   2.138643
    ## 
    ## $Tcm
    ##                  orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGCAGGTGGAT        Tcm       3597         1142   2.224076
    ## AAACCTGGTACGCTGC        Tcm       4175         1185   2.419162
    ## AAACCTGGTCAGTGGA        Tcm       5912         1472   1.962111
    ## AAACCTGGTCGTTGTA        Tcm       8476         1948   2.524776
    ## AAACCTGTCACGAAGG        Tcm       3934         1104   3.050330
    ## AAACCTGTCAGCTGGC        Tcm       4090         1189   3.056235
    ## 
    ## $Tem
    ##                  orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGCAAGGTTCT        Tem       6011         1797  0.9316254
    ## AAACCTGCACGCGAAA        Tem       3740         1186  2.7005348
    ## AAACCTGGTATGGTTC        Tem       1861          776  2.4180548
    ## AAACCTGTCTGGTATG        Tem       2336          955  2.0547945
    ## AAACGGGCATGCCTAA        Tem       6249         1539  2.2083533
    ## AAACGGGGTCTGATTG        Tem       9313         2374  2.2549125
    ## 
    ## $Trm
    ##                  orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGGTCTTGCGG        Trm       4480         1544   2.343750
    ## AAACGGGAGGCTCAGA        Trm       7628         1696   3.224961
    ## AAACGGGGTCAACATC        Trm       4830         1572   2.939959
    ## AAACGGGTCTTGAGAC        Trm       8058         2239   2.593696
    ## AAAGATGGTAGGGTAC        Trm       2378         1003   3.195963
    ## AAAGATGTCCTACAGA        Trm       2934         1027   2.453988

``` r
seulist <- seulist.mt
```

Examine the distribution of nFeature and nCount for each dataset

``` r
seulist %>%
  map_dfr(~.x@meta.data, .id = "seuobj") %>%
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, 
                color = percent.mt))+
    geom_point(size = 1)+
    facet_wrap(~ seuobj, nrow = 2)+
    labs(x = "nCount", y = "nFeature", color = "MT%") +
  theme_classic() +
  geom_hline(yintercept = 200, linetype = 2)
```

![](3_load_transformCounts_files/figure-markdown_github/unnamed-chunk-5-1.png)
Are there any cells with less than 200 expressed genes?

``` r
seulist %>%
  map_dfr(~.x@meta.data, .id = "seuobj") %>%
  ggplot(aes(x = nFeature_RNA))+
    geom_histogram() +
    facet_wrap(~ seuobj, nrow = 2)+
  theme_classic() +
  geom_vline(xintercept = 200, linetype = 2)
```

![](3_load_transformCounts_files/figure-markdown_github/unnamed-chunk-6-1.png)
Summarize as a table

``` r
seulist %>%
  map_dfr(~.x@meta.data, .id = "seuobj") %>%
  group_by(seuobj) %>%
  summarize(min = min(nFeature_RNA),
            mean = mean(nFeature_RNA),
            max = max(nFeature_RNA))
```

    ## # A tibble: 4 × 4
    ##   seuobj   min  mean   max
    ##   <chr>  <int> <dbl> <int>
    ## 1 disTrm   622 1266.  3947
    ## 2 Tcm      686 1259.  2499
    ## 3 Tem      500 1268.  2495
    ## 4 Trm      372 1244.  2442

Are there any cells with more than 5% mitochondrial genes?

``` r
seulist %>%
  map_dfr(~.x@meta.data, .id = "seuobj") %>%
  ggplot(aes(x = percent.mt))+
    geom_histogram() +
    facet_wrap(~ seuobj, nrow = 2)+
  theme_classic() +
  geom_vline(xintercept = 5, linetype = 2)
```

![](3_load_transformCounts_files/figure-markdown_github/unnamed-chunk-8-1.png)
Summarize as a table

``` r
seulist %>%
  map_dfr(~.x@meta.data, .id = "seuobj") %>%
  group_by(seuobj) %>%
  summarize(min = min(percent.mt),
            mean = mean(percent.mt),
            max = max(percent.mt))
```

    ## # A tibble: 4 × 4
    ##   seuobj   min  mean   max
    ##   <chr>  <dbl> <dbl> <dbl>
    ## 1 disTrm 0.826  2.56  4.89
    ## 2 Tcm    0.654  2.44  4.83
    ## 3 Tem    0      2.18  4.95
    ## 4 Trm    0.426  2.35  4.46

## 3. Transform gene counts

*METHOD DETAILS, Single-cell RNA-seq analysis, page e4*: “Gene counts
were scaled to total gene expression and percentage of mitochondrial
genes with a scaling factor of 10,000, and then log-transformed.”

First, take a look at how data are structured inside the filtered seurat
object

``` r
tem <- seulist$Tem

Assays(tem)
```

    ## [1] "RNA"

``` r
tem@assays
```

    ## $RNA
    ## Assay data with 12571 features for 1457 cells
    ## First 10 features:
    ##  Mrpl15, Lypla1, Tcea1, Atp6v1h, Rb1cc1, 4732440D04Rik, Pcmtd1, Gm26901,
    ## Rrs1, Mybl1

There is 1 assay called “RNA”. Inside this, there are 3 slots for
matrices named “counts”, “data”, and “scale.data”. Each matrix consists
of rows = features and columns = cells.

``` r
curr.counts <- tem@assays$RNA@counts
dim(curr.counts) # 12571 features as rows and 1457 cells as columns
```

    ## [1] 12571  1457

``` r
colnames(curr.counts)[1:10]
```

    ##  [1] "AAACCTGCAAGGTTCT" "AAACCTGCACGCGAAA" "AAACCTGGTATGGTTC" "AAACCTGTCTGGTATG"
    ##  [5] "AAACGGGCATGCCTAA" "AAACGGGGTCTGATTG" "AAACGGGTCACCTCGT" "AAACGGGTCTTCAACT"
    ##  [9] "AAAGATGAGTCAATAG" "AAAGATGGTCTTGATG"

``` r
rownames(curr.counts)[1:10]
```

    ##  [1] "Mrpl15"        "Lypla1"        "Tcea1"         "Atp6v1h"      
    ##  [5] "Rb1cc1"        "4732440D04Rik" "Pcmtd1"        "Gm26901"      
    ##  [9] "Rrs1"          "Mybl1"

``` r
curr.data <- tem@assays$RNA@data
dim(curr.data) # 12571 features as rows and 1457 cells as columns
```

    ## [1] 12571  1457

``` r
colnames(curr.data)[1:10]
```

    ##  [1] "AAACCTGCAAGGTTCT" "AAACCTGCACGCGAAA" "AAACCTGGTATGGTTC" "AAACCTGTCTGGTATG"
    ##  [5] "AAACGGGCATGCCTAA" "AAACGGGGTCTGATTG" "AAACGGGTCACCTCGT" "AAACGGGTCTTCAACT"
    ##  [9] "AAAGATGAGTCAATAG" "AAAGATGGTCTTGATG"

``` r
rownames(curr.data)[1:10]
```

    ##  [1] "Mrpl15"        "Lypla1"        "Tcea1"         "Atp6v1h"      
    ##  [5] "Rb1cc1"        "4732440D04Rik" "Pcmtd1"        "Gm26901"      
    ##  [9] "Rrs1"          "Mybl1"

``` r
curr.scale.data <- tem@assays$RNA@scale.data
dim(curr.scale.data) # empty
```

    ## [1] 0 0

Compare how the “counts” and “data” values differ. Pull out values for a
single cell and plot.

``` r
i <- 1
colnames(curr.counts)[i] # cell identifier
```

    ## [1] "AAACCTGCAAGGTTCT"

``` r
df <- data.frame(count = curr.counts[,i], data = curr.data[,i], feature = row.names(curr.counts), row.names = NULL)
head(df)
```

    ##   count data       feature
    ## 1     1    1        Mrpl15
    ## 2     2    2        Lypla1
    ## 3     1    1         Tcea1
    ## 4     1    1       Atp6v1h
    ## 5     1    1        Rb1cc1
    ## 6     0    0 4732440D04Rik

``` r
ggplot(df, aes(x = count, y = data)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey")
```

![](3_load_transformCounts_files/figure-markdown_github/unnamed-chunk-12-1.png)

Normalize the counts and re-examine the data slots. Note that only the
values in the “data” slot have changed.

``` r
tem.s <- NormalizeData(object = tem, normalization.method = "LogNormalize", scale.factor = 10000)

s.counts <- tem.s@assays$RNA@counts
dim(s.counts) # 12571 features as rows and 1457 cells as columns
```

    ## [1] 12571  1457

``` r
s.counts[1:10,1:10]
```

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                  
    ## Mrpl15        1 . . . . 1 . . . .
    ## Lypla1        2 1 . . . 1 . . . .
    ## Tcea1         1 . . 1 1 1 . . . 1
    ## Atp6v1h       1 . . . . . . . . 1
    ## Rb1cc1        1 . . . . . . . . .
    ## 4732440D04Rik . . . 1 . . . . . .
    ## Pcmtd1        . . . . 1 . . . . .
    ## Gm26901       . . . . . . . . . .
    ## Rrs1          . . . . . . . . . 1
    ## Mybl1         . . . . . . . . . .

``` r
#colnames(s.counts)[1:10]
#rownames(s.counts)[1:10]

s.data <- tem.s@assays$RNA@data
dim(s.data) # 12571 features as rows and 1457 cells as columns
```

    ## [1] 12571  1457

``` r
s.data[1:10,1:10]
```

    ## 10 x 10 sparse Matrix of class "dgCMatrix"
    ##                                                                                
    ## Mrpl15        0.9796849 .        . .        .         0.7293672 . . . .        
    ## Lypla1        1.4649284 1.301226 . .        .         0.7293672 . . . .        
    ## Tcea1         0.9796849 .        . 1.664082 0.9556099 0.7293672 . . . 0.8456252
    ## Atp6v1h       0.9796849 .        . .        .         .         . . . 0.8456252
    ## Rb1cc1        0.9796849 .        . .        .         .         . . . .        
    ## 4732440D04Rik .         .        . 1.664082 .         .         . . . .        
    ## Pcmtd1        .         .        . .        0.9556099 .         . . . .        
    ## Gm26901       .         .        . .        .         .         . . . .        
    ## Rrs1          .         .        . .        .         .         . . . 0.8456252
    ## Mybl1         .         .        . .        .         .         . . . .

``` r
#colnames(s.data)[1:10]
#rownames(s.data)[1:10]

s.scale.data <- tem.s@assays$RNA@scale.data
dim(s.scale.data) # empty
```

    ## [1] 0 0

``` r
## plot
i <- 1
colnames(s.counts)[i] # cell identifier
```

    ## [1] "AAACCTGCAAGGTTCT"

``` r
df <- data.frame(count = s.counts[,i], data = s.data[,i], feature = row.names(s.counts), row.names = NULL)
ggplot(df, aes(x = count, y = data)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey")
```

![](3_load_transformCounts_files/figure-markdown_github/unnamed-chunk-13-1.png)

Normalize each of the seurat objects in the list.

Here is an example of how to do this with a for-loop

``` r
sampleNames <- names(seulist)
seulist.n <- list()
for(i in 1:length(sampleNames)){
  seulist.n[[i]] <- NormalizeData(object = seulist[[i]], normalization.method = "LogNormalize", 
                                  scale.factor = 10000)
}
names(seulist.n) <- sampleNames
```

Here is an example of how to this using purrr::map and pipes

``` r
# seulist %>%
#   map(~NormalizeData(object = .x, normalization.method = "LogNormalize", scale.factor = 10000)) -> seulist.n
```
