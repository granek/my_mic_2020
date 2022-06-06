-   [Overview](#overview)
-   [Load data](#load-data)
-   [Highly variable gene (HVG)
    identification](#highly-variable-gene-hvg-identification)
-   [Dimension reduction](#dimension-reduction)
-   [Clustering and visualization](#clustering-and-visualization)
    -   [Apply to the list of seurat
        objects](#apply-to-the-list-of-seurat-objects)
    -   [Examine the tSNE plots
        (Fig2A).](#examine-the-tsne-plots-fig2a.)
-   [Cell-type specific marker gene
    detection](#cell-type-specific-marker-gene-detection)

Before we start, we set up the R environment.

``` r
knitr::opts_chunk$set(echo = T, message = FALSE, warning = FALSE)
stdt <- date()

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

# Input paths
#wd <- "~/Documents/Teaching/chsi-mic-2022"
wd <- "/hpc/group/chsi-mic-2022"
intermeddir <- file.path(wd, "intermed")

# Output paths
scratchdir <- "~/Documents/Temporary/scratch"
# scratchdir <- "/work"
# username <- Sys.info()[["user"]]
# outdir <- file.path(scratchdir,username,"output")
# dir.create(outdir) # uncomment if not already there
```

## Overview

**Data source:** Christian et al. 2021 Cell Reports
(<https://doi.org/10.1016/j.celrep.2021.109118>)

**Goal of this workshop:** Illustrate the clustering and cell-type
annotation procedures.

Note that the Seurat - [Guided Clustering
Tutorial](https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
includes a lot of overlapping information.

**METHOD DETAILS, Single-cell RNA-seq analysis, page e4** -
Pre-processing: “Raw short reads were demultiplexed, filtering and
mapped to mouse genome GRCm38/mm10 using cellranger v2.02.” - QC: “The
gene count matrices from cellranger were subjected to quality control,
pre-processing and clustering using the R Seurat 2.3.4 package (Butler
et al., 2018). Low-quality cells that had less than 200 expressed genes
and more than 5% mitochondrial genes were filtered out.” - Analysis:
“Gene counts were scaled to total gene expression and percentage of
mitochondrial genes with a scaling factor of 10,000, and then
log-transformed. The high dimensional data for each sample were reduced
by PCA and t-Distributed Stochastics Neighbor Embedding (tSNE). We used
the FindCluster function to group clusters in each sample with a
resolution of 0.6. Differential expressed genes (DEGs) were identified
using the Wilcoxon rank-sum test.”

**What’s covered in this workshop:** - Highly variable gene
identification - Dimension reduction - Clustering and visualization -
Cell-type specific marker gene detection - How to generate the figures
for the paper?

## Load data

First, we load the normalized UMI count single-cell data.

``` r
seulist <- readRDS(file.path(intermeddir, "seulist_drc.rds"))
ls()
```

    ## [1] "intermeddir" "scratchdir"  "seulist"     "stdt"        "wd"

``` r
seulist
```

    ## $disTrm
    ## An object of class Seurat 
    ## 12064 features across 863 samples within 1 assay 
    ## Active assay: RNA (12064 features, 0 variable features)
    ##  2 dimensional reductions calculated: pca, tsne
    ## 
    ## $Tcm
    ## An object of class Seurat 
    ## 12840 features across 4009 samples within 1 assay 
    ## Active assay: RNA (12840 features, 0 variable features)
    ##  2 dimensional reductions calculated: pca, tsne
    ## 
    ## $Tem
    ## An object of class Seurat 
    ## 12571 features across 1457 samples within 1 assay 
    ## Active assay: RNA (12571 features, 0 variable features)
    ##  2 dimensional reductions calculated: pca, tsne
    ## 
    ## $Trm
    ## An object of class Seurat 
    ## 11778 features across 459 samples within 1 assay 
    ## Active assay: RNA (11778 features, 0 variable features)
    ##  2 dimensional reductions calculated: pca, tsne

``` r
tem <- seulist[["Tem"]]
```

## Highly variable gene (HVG) identification

**First, find the top 2000 most variable gene features.** This function
finds the features, stashes them inside the seurat object, then returns
the seurat object.

``` r
tem.vf <- FindVariableFeatures(tem, selection.method = "vst", nfeatures = 2000)
```

You can extract the variable features from the seurat object using these
methods.

``` r
VariableFeatures(tem.vf)[1:10] # option 1
```

    ##  [1] "Hist1h2ap" "Gzma"      "Ccl4"      "Ccl3"      "Ifitm1"    "Ifng"     
    ##  [7] "Hmgb2"     "Actb"      "Ccl5"      "Hist1h2ae"

``` r
tem.vf@assays$RNA@var.features[1:10] # option 2
```

    ##  [1] "Hist1h2ap" "Gzma"      "Ccl4"      "Ccl3"      "Ifitm1"    "Ifng"     
    ##  [7] "Hmgb2"     "Actb"      "Ccl5"      "Hist1h2ae"

``` r
# note that this slot is empty in the seurat object prior to the FindVariableFeatures() command
tem@assays$RNA@var.features
```

    ## logical(0)

In the `FindVariableFeatures()` function, we used `vst` as the selection
method. Intermediate data associated with that selection method is
stored inside the seurat object here:

``` r
head(tem.vf@assays$RNA@meta.features)
```

    ##                 vst.mean vst.variance vst.variance.expected
    ## Mrpl15        0.24021963   0.25818896            0.28347872
    ## Lypla1        0.17638984   0.20993668            0.20285509
    ## Tcea1         0.37474262   0.48721877            0.46985696
    ## Atp6v1h       0.15168154   0.17271867            0.17290937
    ## Rb1cc1        0.13452299   0.14672630            0.15227544
    ## 4732440D04Rik 0.01990391   0.01952114            0.02134118
    ##               vst.variance.standardized vst.variable
    ## Mrpl15                        0.9107878        FALSE
    ## Lypla1                        1.0349096        FALSE
    ## Tcea1                         1.0369513        FALSE
    ## Atp6v1h                       0.9988971        FALSE
    ## Rb1cc1                        0.9635585        FALSE
    ## 4732440D04Rik                 0.9147173        FALSE

Seurat has a function to generate a plot to show how variable these
“variable features” are.

``` r
top10 <- head(VariableFeatures(tem.vf), 10)
plot1 <- VariableFeaturePlot(tem.vf)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
plot2
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-4-2.png)

## Dimension reduction

*METHOD DETAILS, Single-cell RNA-seq analysis, page e4*: “The high
dimensional data for each sample were reduced by PCA…”

Specifically, shift the expression of each gene, so that the mean
expression across cells is 0 and the variance across cells is 1. If we
didn’t do this, highly-expressed and highly-variable genes would
dominate the outcome of down-stream analyses.

``` r
all.genes <- rownames(tem.vf)
tem.sAll <- ScaleData(tem.vf, features = all.genes) 
tem.sVF <- ScaleData(tem.vf) 

# note that if you do not specify features as all.genes, the default behavior is to only use the 2000 variable features downstream
```

Examine the scaled data.

``` r
# prior to scaling
dim(tem.vf@assays$RNA@scale.data) # empty
```

    ## [1] 12571  1457

``` r
# after scaling with all genes
dim(tem.sAll@assays$RNA@scale.data) # 12571 features as rows and 1457 cells as columns
```

    ## [1] 12571  1457

``` r
# after scaling with the 2000 variable features
dim(tem.sVF@assays$RNA@scale.data) # 2000 features as rows and 1457 cells as columns
```

    ## [1] 2000 1457

We perform PCA on the scaled data.

``` r
tem.pAll <- RunPCA(tem.sAll, features = all.genes) # this takes a minute or so...
tem.pVF <- RunPCA(tem.sVF)
# again, note that the default behavior is to use only the variables to compute the PCA
```

The `RunPCA()` function automatically spits out information about how
features loaded into PC space. But like most Seurat functions, you will
need to look inside the object to get details.

Here is where you will find detailed info about the reduction you
performed

``` r
tem.pVF@reductions # the seurat object keeps track of reductions you have performed here
```

    ## $pca
    ## A dimensional reduction object with key PC_ 
    ##  Number of dimensions: 50 
    ##  Projected dimensional reduction calculated:  FALSE 
    ##  Jackstraw run: FALSE 
    ##  Computed using assay: RNA 
    ## 
    ## $tsne
    ## A dimensional reduction object with key tSNE_ 
    ##  Number of dimensions: 2 
    ##  Projected dimensional reduction calculated:  FALSE 
    ##  Jackstraw run: FALSE 
    ##  Computed using assay: RNA

``` r
names(tem.pVF@commands) # more broadly, there is a log of all commands stored here
```

    ## [1] "NormalizeData.RNA"        "FindNeighbors.RNA.pca"   
    ## [3] "FindClusters"             "RunTSNE"                 
    ## [5] "FindVariableFeatures.RNA" "ScaleData.RNA"           
    ## [7] "RunPCA.RNA"

``` r
tem.pVF@commands$RunPCA.RNA # with even more details
```

    ## Command: RunPCA(tem.sVF)
    ## Time: 2022-06-05 21:12:45
    ## assay : RNA 
    ## npcs : 50 
    ## rev.pca : FALSE 
    ## weight.by.var : TRUE 
    ## verbose : TRUE 
    ## ndims.print : 1 2 3 4 5 
    ## nfeatures.print : 30 
    ## reduction.name : pca 
    ## reduction.key : PC_ 
    ## seed.use : 42

Next, we can use the visualization tools to the choose the proper number
of PCs.

``` r
VizDimLoadings(tem.pVF, dims = 1:2, reduction = "pca")
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
DimPlot(tem.pVF, reduction = "pca")
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-9-2.png)

``` r
DimHeatmap(tem.pVF, dims = 1, cells = 500, balanced = TRUE)
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-9-3.png)

Check the first 15 PCs

``` r
DimHeatmap(tem.pVF, dims = 1:15, cells = 500, balanced = TRUE)
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-10-1.png)

To refine the number of dimensions, we can use…

**Option 1:** A resampling test inspired by the JackStraw procedure
developed by the Seurat authors. The idea is..”to permute a subset of
the data (1% by default) and rerun PCA, constructing a ‘null
distribution’ of feature scores, and repeat this procedure. We identify
‘significant’ PCs as those who have a strong enrichment of low p-value
features….The `JackStrawPlot()` function provides a visualization tool
for comparing the distribution of p-values for each PC with a uniform
distribution (dashed line). ‘Significant’ PCs will show a strong
enrichment of features with low p-values (solid curve above the dashed
line).”(<https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>).

``` r
#tem.js <- JackStraw(tem.pVF, num.replicate = 100) # this takes a minute or two, so we can just load a version that we previously created, cooking-show style.

# Students will not be able to run the saveRDS command because 
# everything the data/ and intermed/ folders will need to be write-protected
#saveRDS(tem.js, file = file.path(intermeddir, "tem_js.rds"))
tem.js <- readRDS(file.path(intermeddir, "tem_js.rds"))

tem.js <- ScoreJackStraw(tem.js, dims = 1:20)
JackStrawPlot(tem.js, dims = 1:20)
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-11-1.png)

There seems to be a drop-off in significance somewhere between PC5 and
PC10 but it is difficult to tell with the colors. This suggests we
should use \~10 dimensions.

**Option 2:** We can also just take a look at the percentage of variance
explained by each additional PC. We don’t want to bother including PCs
that explain very little additional variance so we can look for the
“elbow” where the variance explained (or inversely, the standard
deviation) levels-off. Seurat’s tool, `ElbowPlot()` displays standard
deviation on the y-axis.

``` r
ElbowPlot(tem.pVF)
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-12-1.png)
Here we can see an “elbow” that forms around PC10 to PC15, suggesting
that we should to use at least that many dimensions and, beyond that, we
start seeing diminishing returns.

## Clustering and visualization

*METHOD DETAILS, Single-cell RNA-seq analysis, page 4e*: ” ..and
t-Distributed Stochastics Neighbor Embedding (tSNE). We used the
`FindCluster()` function to group clusters in each sample with a
resolution of 0.6.”

``` r
tem.pVF <- FindNeighbors(tem.pVF, dims = 1:10)
tem.pVF <- FindClusters(tem.pVF, resolution = 0.6)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 1457
    ## Number of edges: 50426
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8143
    ## Number of communities: 6
    ## Elapsed time: 0 seconds

The cluster IDs are stored here.

``` r
Idents(tem.pVF)[1:10]
```

    ## AAACCTGCAAGGTTCT AAACCTGCACGCGAAA AAACCTGGTATGGTTC AAACCTGTCTGGTATG 
    ##                5                0                3                3 
    ## AAACGGGCATGCCTAA AAACGGGGTCTGATTG AAACGGGTCACCTCGT AAACGGGTCTTCAACT 
    ##                0                4                1                0 
    ## AAAGATGAGTCAATAG AAAGATGGTCTTGATG 
    ##                4                4 
    ## Levels: 0 1 2 3 4 5

After finding clusters, project into 2D for visualization using tSNE.

``` r
tem.pVf.tsne <- RunTSNE(tem.pVF, dims = 1:10) 
# note that the seurat authors suggest using the same PCs as input to the clustering analysis
DimPlot(tem.pVf.tsne, reduction = "tsne") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-15-1.png)

Another option is to do this with UMAP instead of tSNE. They should look
pretty similar.

``` r
tem.pVf.umap <- RunUMAP(tem.pVF, dims = 1:10)
DimPlot(tem.pVf.umap, reduction = "umap") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-16-1.png)

### Apply to the list of seurat objects

Make a function that contains that takes a seurat object as input, does
all the operations we want, then returned the updated seurat object.

``` r
dimRed_and_cluster <- function(seuobj){
  
  #seuobj <- FindVariableFeatures(seuobj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seuobj)
  seuobj <- ScaleData(seuobj, features = all.genes) 
  seuobj <- RunPCA(seuobj, features = all.genes)
  seuobj <- FindNeighbors(seuobj, dims = 1:10)
  seuobj <- FindClusters(seuobj, resolution = 0.6)
  seuobj <- RunTSNE(seuobj, dims = 1:10)
  
  return(seuobj)

}
```

Apply that function to each seurat object in the list of objects that we
normalized

``` r
#seulist.n %>%
#  map(~dimRed_and_cluster(.x)) -> seulist.drc
# this takes a minute or two, so we can just load a version that we previously created, cooking-show style.

# Students will not be able to run the saveRDS command because 
# everything the data/ and intermed/ folders will need to be write-protected
#saveRDS(seulist.drc, file = file.path(intermeddir, "seulist_drc.rds"))
seulist.drc <- readRDS(file.path(intermeddir, "seulist_drc.rds"))
```

### Examine the tSNE plots (Fig2A).

These should resemble plots in Fig2, page 5.

**Fig 2A**. They are similar but have different orientations and numbers
of clusters.

``` r
DimPlot(seulist.drc[["Tem"]], reduction = "tsne") + ggtitle("Tem") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
# Clusters 0, 1, 3 correspond to p1?
# Cluster 3 corresponds to p2?
# Cluster 4 corresponds to p3?
# Cluster 5 corresponds to p4?

DimPlot(seulist.drc[["Tcm"]], reduction = "tsne") + ggtitle("Tcm") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-19-2.png)

``` r
DimPlot(seulist.drc[["Trm"]], reduction = "tsne") + ggtitle("Trm") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-19-3.png)

``` r
DimPlot(seulist.drc[["disTrm"]], reduction = "tsne") + ggtitle("disTrm") + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-19-4.png)

Overlay the expression of 1 gene on top of the Tem cells in tSNE space
as shown in **Fig 2B**.

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Ifng")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
#FeaturePlot(seulist.drc[["Tem"]], features = c("GzmB"))
#The following requested variables were not found: GzmB????
FeaturePlot(seulist.drc[["Tem"]], features = c("Klrc1")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-20-2.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Nkg7")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-20-3.png)

Overlay … as shown in **Fig 2C**.

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Runx3")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Id2")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-21-2.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Id3")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-21-3.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Bcl2")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-21-4.png)

Overlay … as shown in **Fig 2D**.

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Ccnb2")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Cdk1")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-22-2.png)

``` r
FeaturePlot(seulist.drc[["Tem"]], features = c("Mki67")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-22-3.png)

## Cell-type specific marker gene detection

*METHOD DETAILS, Single-cell RNA-seq analysis, page 4e*: “Differential
expressed genes (DEGs) were identified using the Wilcoxon rank-sum
test.”

#### Identify Tem cluster biomarkers (Fig2B-D)

``` r
tem <- seulist.drc[["Tem"]]
tem@meta.data %>%
  pull(seurat_clusters) %>% unique()
```

    ## [1] 5 1 0 4 3 2
    ## Levels: 0 1 2 3 4 5

For the purpose of trying to reproduce the manuscript results, we will
assume that the authors use defaults. Here are a few default arguments
to be aware of:

**min.pct = 0.1** : Do not test genes that make up fewer than 0.1
fraction (10% ?) of the total reads in either of the populations tested.
Meant to speed up the function by not testing genes that are very
infrequently expressed. Default is 0.1 – **NOTE: Marissa tried to
re-write this so that it is clearer, please check**

**max.cells.per.ident = Inf** : “This will downsample each identity
class to have no more cells than whatever this is set to. While there is
generally going to be a loss in power, the speed increases can be
significant and the most highly differentially expressed features will
likely still rise to the top”
(<https://satijalab.org/seurat/articles/pbmc3k_tutorial.html>)

**logfc.threshold = 0.25** “Limit testing to genes which show, on
average, at least X-fold difference (log-scale) between the two groups
of cells. Default is 0.25 Increasing logfc.threshold speeds up the
function, but can miss weaker signals.” (Seurat help)

**test.use = “wilcox** Denotes which test to use; see Seurat help for
alternatives.

**min.cells.group = 3** “Minimum number of cells in one of the groups”
(Seurat help)

**ident.2** “A second identity class for comparison; if NULL, use all
other cells for comparison” (Seurat help)

Cluster 0

``` r
cluster0.markers <- FindMarkers(tem, ident.1 = 0)
cluster0.markers %>%
  arrange(-avg_log2FC) %>%
  head(n = 5)
```

    ##                p_val avg_log2FC pct.1 pct.2    p_val_adj
    ## Bcl2    2.175929e-28  0.8085540 0.855 0.608 2.735360e-24
    ## Dapl1   3.665931e-09  0.6328763 0.126 0.039 4.608442e-05
    ## Rps28   3.679618e-70  0.5631877 1.000 0.999 4.625647e-66
    ## Rps27rt 8.884549e-46  0.5063820 1.000 0.996 1.116877e-41
    ## Gm10260 3.618803e-44  0.5046848 1.000 0.994 4.549197e-40

Check the expression of the feature genes.

``` r
FeaturePlot(tem, features = c("Bcl2")) + coord_equal()
```

![](4_clustering_annotation_files/figure-markdown_github/unnamed-chunk-25-1.png)
