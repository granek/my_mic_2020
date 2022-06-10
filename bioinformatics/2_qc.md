-   [Convert cellranger output to a Seurat
    object](#convert-cellranger-output-to-a-seurat-object)
-   [Calculate mitrochondiral gene
    percentage](#calculate-mitrochondiral-gene-percentage)
-   [Learn data quality attributes, Initial
    QC](#learn-data-quality-attributes-initial-qc)
-   [Exclude data from low-quality
    cells](#exclude-data-from-low-quality-cells)
    -   [nCount](#ncount)
    -   [nFeature](#nfeature)
    -   [MT%](#mt)
    -   [Do the filtering](#do-the-filtering)
-   [Final QC](#final-qc)
-   [Session Info](#session-info)

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
# Paths
cr.outs <- "/hpc/group/chsi-mic-2022/data/pbmc8k/pbmc8k/outs"
intermeddir <- "/hpc/group/chsi-mic-2022/intermed"
```

**Goal of this workshop:** Learn how to QC scRNAseq data

**What’s covered in this workshop:** - Convert cellranger output to a
Seurat object - Calculate mitrochondiral gene percentage - Learn data
quality attributes - Exclude data from low-quality cells

**Data source:** scRNA-seq dataset of 8K human peripheral blood
mononuclear cells (PBMCs) freely available from 10X Genomics
(<https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc8k>)

## Convert cellranger output to a Seurat object

We can use code like this to load the cell ranger output and store it in
a Seurat object.

``` r
Read10X_h5(file.path(cr.outs, "filtered_feature_bc_matrix.h5")) %>% 
  CreateSeuratObject() -> pbmc8k
pbmc8k
```

    ## An object of class Seurat 
    ## 36601 features across 8788 samples within 1 assay 
    ## Active assay: RNA (36601 features, 0 variable features)

Alternatively, you can load all the files inside the
`filtered_feature_bc_matrix` folder Here are the contents of this folder

``` r
dir(file.path(cr.outs, "filtered_feature_bc_matrix"))
```

    ## [1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"

Here is how you can load the data into a Seurat object.

``` r
Read10X(file.path(cr.outs, "filtered_feature_bc_matrix")) %>% 
  CreateSeuratObject() -> pbmc8k.alt
pbmc8k.alt
```

    ## An object of class Seurat 
    ## 36601 features across 8788 samples within 1 assay 
    ## Active assay: RNA (36601 features, 0 variable features)

Since these objects are large, I am going to remove the alternative
object from the R environment

``` r
rm(pbmc8k.alt)
```

## Calculate mitrochondiral gene percentage

First, take look at the existing cell-level metadata in one of the
seurat objects

``` r
head(pbmc8k@meta.data)
```

    ##                       orig.ident nCount_RNA nFeature_RNA
    ## AAACCTGAGCATCATC-1 SeuratProject       3619         1631
    ## AAACCTGAGCTAACTC-1 SeuratProject       2200         1169
    ## AAACCTGAGCTAGTGG-1 SeuratProject       5472         1917
    ## AAACCTGCACATTAGC-1 SeuratProject       3374         1343
    ## AAACCTGCACTGTTAG-1 SeuratProject       5941         2241
    ## AAACCTGCATAGTAAG-1 SeuratProject       5369         2048

orig.ident = sample name nCount_RNA = number of reads per cell
nFeature_RNA = number of gene features per cell

The rownames are gene symbols. Investigating the names a bit more, the
mitochondrial genes can be identified by the prefix “mt-”

``` r
curr.count <- pbmc8k@assays$RNA@counts
rownames(curr.count)[1:10]
```

    ##  [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3" 
    ##  [6] "AL627309.2"  "AL627309.5"  "AL627309.4"  "AP006222.2"  "AL732372.1"

``` r
rownames(curr.count)[grepl("^MT-", rownames(curr.count))] # these look like mitochondrial genes!
```

    ##  [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3" 
    ##  [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

Add the mitochondrial percentage for 1 seurat object

``` r
head(pbmc8k@meta.data)
```

    ##                       orig.ident nCount_RNA nFeature_RNA
    ## AAACCTGAGCATCATC-1 SeuratProject       3619         1631
    ## AAACCTGAGCTAACTC-1 SeuratProject       2200         1169
    ## AAACCTGAGCTAGTGG-1 SeuratProject       5472         1917
    ## AAACCTGCACATTAGC-1 SeuratProject       3374         1343
    ## AAACCTGCACTGTTAG-1 SeuratProject       5941         2241
    ## AAACCTGCATAGTAAG-1 SeuratProject       5369         2048

``` r
pbmc8k[["percent.mt"]] <- PercentageFeatureSet(pbmc8k, pattern = "^MT-")
head(pbmc8k@meta.data)
```

    ##                       orig.ident nCount_RNA nFeature_RNA percent.mt
    ## AAACCTGAGCATCATC-1 SeuratProject       3619         1631   2.625035
    ## AAACCTGAGCTAACTC-1 SeuratProject       2200         1169   4.363636
    ## AAACCTGAGCTAGTGG-1 SeuratProject       5472         1917   1.608187
    ## AAACCTGCACATTAGC-1 SeuratProject       3374         1343   1.215175
    ## AAACCTGCACTGTTAG-1 SeuratProject       5941         2241   2.794142
    ## AAACCTGCATAGTAAG-1 SeuratProject       5369         2048   3.054573

## Learn data quality attributes, Initial QC

Seurat has a built-in function to plot the three key aspects of data
quality that we need to check. Below we will explain each attribute and
show how to generate similar plots using ggplot.

``` r
VlnPlot(pbmc8k, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-8-1.png)

Number of molecules aka reads detected per cell (nCount)

-   A very low number of reads per cell could indicate a sequencing
    failure.
-   A very high number of reads per cell could indicate more than one
    cell was actually sequenced.

``` r
ggplot(pbmc8k@meta.data, aes(x = orig.ident, y = nCount_RNA))+
    geom_jitter(shape = 16, position = position_jitter(0.2))+
    geom_violin(trim = F, alpha = 0.7) +
  labs(x = "Group", y = "nCount") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        strip.text = element_text(size = 6))
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-9-1.png)

Number of gene features detected per cell (nFeature)

-   A very low number of gene features per cell could indicate a library
    prep or sequencing failure.

``` r
ggplot(pbmc8k@meta.data, aes(x = orig.ident, y = nFeature_RNA))+
    geom_jitter(shape = 16, position = position_jitter(0.2))+
    geom_violin(trim = F, alpha = 0.7)+
  labs(x = "Group", y = "nFeature") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        strip.text = element_text(size = 6))
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-10-1.png)

Percentage of mitochondrial genes per cell (MT%)

-   A high percentage of mitochondrial genes (MT%) indicates a cell may
    be dead or dying based on the expectation that, if a cell is
    ruptured, non-MT genes will leak out first and increase the relative
    abundance of MT genes.

``` r
ggplot(pbmc8k@meta.data, aes(x = orig.ident, y = percent.mt)) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
    geom_violin(trim=FALSE, alpha=0.7)+
  scale_y_continuous(limits = c(0, 10)) +
  labs(x = "Group", y = "MT%") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust=1),
        strip.text = element_text(size = 6)) +
  geom_hline(yintercept = 5, linetype = 2)
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-11-1.png)

Correlation between nCounts and nFeatures, colored by MT%

-   The number of gene features detected in a cell (nFeatures) tends to
    increase with library size (nCounts). Divergence from this
    correlation could indicate low-quality data, e.g. often observed in
    high MT% cells.

``` r
ggplot(pbmc8k@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, 
                color = percent.mt))+
    geom_point(size = 1)+
    labs(x = "nCount", y = "nFeature", color = "MT%") +
  theme_classic() +
  geom_hline(yintercept = 200, linetype = 2)
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-12-1.png)

## Exclude data from low-quality cells

“Low-quality cells that had less than 200 expressed genes and more than
5% mitochondrial genes were filtered out.”

### nCount

``` r
cl <- 6 # min number of total reads per cell (nCount is analogous to library size or total number of reads)
ch <- 11 # max number of total reads per cell
```

Filter out cells with log(nCount) \> `ch` & \< `cl`

-   Note: Different y-axis scales

``` r
ggplot(pbmc8k@meta.data, aes(x = log(nCount_RNA))) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = c(cl, ch), 
             color = "blue", linetype = "dashed") +
  labs(x = "log(nCount)", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  theme_classic()
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-14-1.png)

### nFeature

``` r
fl <- 5.5 # min number of gene features per cell
```

Filtering out cells with log(nFeature) \< `fl`

-   Note: Different y-axis scales

``` r
ggplot(pbmc8k@meta.data, aes(x = log(nFeature_RNA))) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = fl, 
             color = "blue", linetype = "dashed") +
  scale_x_continuous(breaks = c(seq(0, 20, 1), fl)) +
  labs(x = "log(nCount)", y = "Frequency") +
  theme_classic()
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-16-1.png)

### MT%

``` r
data.frame(threshold = 1:40) %>%
  full_join(pbmc8k@meta.data, by = character()) %>%
  mutate(type = ifelse(percent.mt < threshold, "Kept", 
                       ifelse(percent.mt >= threshold, "Removed", NA))) %>%
  group_by(orig.ident, threshold, type) %>%
  summarize(ncells = n()) -> tmp

ml <- 5
ggplot(tmp, aes(x=threshold, y=ncells, colour=type))+
    geom_line() + 
  geom_vline(xintercept = ml, color = "blue", linetype = "dashed") + 
  scale_x_continuous(breaks = seq(0, 40, by = 5)) +
    labs(x="Thresholds for MT%", y="Number of Cells")+
  facet_wrap(~orig.ident, nrow = 3, 
             scales = "free_y")+
  theme_classic()
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-17-1.png)

Filtering out cells with MT% \> `ml` There aren’t any of these cells
because they have already been filtered.

``` r
ggplot(pbmc8k@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = ml,
             color = "blue", linetype = "dashed") +
  # scale_x_continuous(breaks = c(seq(2, 14, 2), fl)) +
  labs(x = "MT%", y = "Frequency") +
  theme_classic()
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-18-1.png)

### Do the filtering

Pre

``` r
nrow(pbmc8k@meta.data)
```

    ## [1] 8788

Filter

``` r
pbmc8k.filt <- subset(pbmc8k, nCount_RNA > exp(cl) & nCount_RNA < exp(ch) &
                        nFeature_RNA > exp(fl) & 
                        percent.mt < ml)
```

Post

``` r
nrow(pbmc8k.filt@meta.data)
```

    ## [1] 8456

Number of cells removed by filtering

``` r
nrow(pbmc8k@meta.data) - nrow(pbmc8k.filt@meta.data)
```

    ## [1] 332

## Final QC

nCount

``` r
ggplot(pbmc8k.filt@meta.data, aes(x = orig.ident, y = nCount_RNA))+
    geom_jitter(shape = 16, position = position_jitter(0.2))+
    geom_violin(trim = F, alpha = 0.7) +
  labs(x = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-23-1.png)

nFeature

``` r
ggplot(pbmc8k.filt@meta.data, aes(x = orig.ident, y = nFeature_RNA)) +
    geom_jitter(shape = 16, position = position_jitter(0.2))+
    geom_violin(trim = F, alpha = 0.7)+
  labs(x = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust=1))
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-24-1.png)

MT%

``` r
ggplot(pbmc8k.filt@meta.data, aes(x = orig.ident, y = percent.mt)) +
  geom_jitter(shape=16, position=position_jitter(0.2))+
    geom_violin(trim=FALSE, alpha=0.7)+
  labs(x = "Group", y = "MT%") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-25-1.png)

Correlation between nCounts and nFeatures, colored by MT%

``` r
ggplot(pbmc8k.filt@meta.data, aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt))+
    geom_point(size = 1)+
    labs(x = "nCount", y = "nFeature", color = "MT%") +
  theme_classic()
```

![](2_qc_files/figure-markdown_github/unnamed-chunk-26-1.png)

Save the filtered seurat object as a R data object

``` r
newfile <- file.path(intermeddir, "pbmc8k-filt.rds")
#saveRDS(pbmc8k.filt, newfile)

tools::md5sum(newfile)
```

    ## /hpc/group/chsi-mic-2022/intermed/pbmc8k-filt.rds 
    ##                "d12aac6e3486482432acc3dc33121197"

## Session Info

``` r
sessionInfo()
```

    ## R version 4.2.0 (2022-04-22)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 20.04.4 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] sp_1.4-7           SeuratObject_4.1.0 Seurat_4.1.1.9001  forcats_0.5.1     
    ##  [5] stringr_1.4.0      dplyr_1.0.9        purrr_0.3.4        readr_2.1.2       
    ##  [9] tidyr_1.2.0        tibble_3.1.7       ggplot2_3.3.6      tidyverse_1.3.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] readxl_1.4.0          backports_1.4.1       plyr_1.8.7           
    ##   [4] igraph_1.3.1          lazyeval_0.2.2        splines_4.2.0        
    ##   [7] listenv_0.8.0         scattermore_0.8       digest_0.6.29        
    ##  [10] htmltools_0.5.2       fansi_1.0.3           magrittr_2.0.3       
    ##  [13] tensor_1.5            cluster_2.1.3         ROCR_1.0-11          
    ##  [16] tzdb_0.3.0            globals_0.15.0        modelr_0.1.8         
    ##  [19] matrixStats_0.62.0    spatstat.sparse_2.1-1 colorspace_2.0-3     
    ##  [22] rvest_1.0.2           ggrepel_0.9.1         haven_2.5.0          
    ##  [25] xfun_0.31             crayon_1.5.1          jsonlite_1.8.0       
    ##  [28] progressr_0.10.1      spatstat.data_2.2-0   survival_3.2-13      
    ##  [31] zoo_1.8-10            glue_1.6.2            polyclip_1.10-0      
    ##  [34] gtable_0.3.0          leiden_0.4.2          future.apply_1.9.0   
    ##  [37] abind_1.4-5           scales_1.2.0          DBI_1.1.2            
    ##  [40] spatstat.random_2.2-0 miniUI_0.1.1.1        Rcpp_1.0.8.3         
    ##  [43] viridisLite_0.4.0     xtable_1.8-4          reticulate_1.25      
    ##  [46] spatstat.core_2.4-4   bit_4.0.4             htmlwidgets_1.5.4    
    ##  [49] httr_1.4.3            RColorBrewer_1.1-3    ellipsis_0.3.2       
    ##  [52] ica_1.0-2             farver_2.1.0          pkgconfig_2.0.3      
    ##  [55] uwot_0.1.11           dbplyr_2.1.1          deldir_1.0-6         
    ##  [58] utf8_1.2.2            labeling_0.4.2        tidyselect_1.1.2     
    ##  [61] rlang_1.0.2           reshape2_1.4.4        later_1.3.0          
    ##  [64] munsell_0.5.0         cellranger_1.1.0      tools_4.2.0          
    ##  [67] cli_3.3.0             generics_0.1.2        broom_0.8.0          
    ##  [70] ggridges_0.5.3        evaluate_0.15         fastmap_1.1.0        
    ##  [73] yaml_2.3.5            goftest_1.2-3         knitr_1.39           
    ##  [76] bit64_4.0.5           fs_1.5.2              fitdistrplus_1.1-8   
    ##  [79] RANN_2.6.1            pbapply_1.5-0         future_1.26.1        
    ##  [82] nlme_3.1-157          mime_0.12             ggrastr_1.0.1        
    ##  [85] xml2_1.3.3            hdf5r_1.3.5           compiler_4.2.0       
    ##  [88] rstudioapi_0.13       beeswarm_0.4.0        plotly_4.10.0        
    ##  [91] png_0.1-7             spatstat.utils_2.3-1  reprex_2.0.1         
    ##  [94] stringi_1.7.6         highr_0.9             rgeos_0.5-9          
    ##  [97] lattice_0.20-45       Matrix_1.4-1          vctrs_0.4.1          
    ## [100] pillar_1.7.0          lifecycle_1.0.1       spatstat.geom_2.4-0  
    ## [103] lmtest_0.9-40         RcppAnnoy_0.0.19      data.table_1.14.2    
    ## [106] cowplot_1.1.1         irlba_2.3.5           httpuv_1.6.5         
    ## [109] patchwork_1.1.1       R6_2.5.1              promises_1.2.0.1     
    ## [112] KernSmooth_2.23-20    gridExtra_2.3         vipor_0.4.5          
    ## [115] parallelly_1.31.1     codetools_0.2-18      MASS_7.3-57          
    ## [118] assertthat_0.2.1      withr_2.5.0           sctransform_0.3.3    
    ## [121] mgcv_1.8-40           parallel_4.2.0        hms_1.1.1            
    ## [124] grid_4.2.0            rpart_4.1.16          rmarkdown_2.14       
    ## [127] Rtsne_0.16            shiny_1.7.1           lubridate_1.8.0      
    ## [130] ggbeeswarm_0.6.0

``` r
print(paste("Start Time:  ",stdt))
```

    ## [1] "Start Time:   Fri Jun 10 17:45:03 2022"

``` r
print(paste("End Time:  ",date()))
```

    ## [1] "End Time:   Fri Jun 10 17:45:30 2022"
