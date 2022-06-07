-   [Load data](#load-data)
-   [1. Merge Seurat objects](#merge-seurat-objects)
-   [2. Create a Monocle object](#create-a-monocle-object)
-   [3. Run tranjectory analysis](#run-tranjectory-analysis)
-   [4. Test for DE over psuedotime](#test-for-de-over-psuedotime)
-   [5. Additional resources](#additional-resources)

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
library(monocle); package.version("monocle")
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: VGAM

    ## Loading required package: stats4

    ## Loading required package: splines

    ## 
    ## Attaching package: 'VGAM'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     fill

    ## Loading required package: DDRTree

    ## Loading required package: irlba

    ## [1] "2.24.0"

``` r
# Paths
wd <- "/hpc/group/chsi-mic-2022"
intermeddir <- file.path(wd, "intermed")
outdir <- file.path("/work",Sys.info()[["user"]],"output")
#outdir
```

**Goal of this workshop:** Learn to analyze cell trajectories with
monocle

**What’s covered in this workshop:**

-   Merge Seurat objects
-   Create a Monocle object
-   Run trajectory analysis
-   Test for DE over pseudotime
-   Additional resources

**Data source:** Christian et al. 2021 Cell Reports
(<https://doi.org/10.1016/j.celrep.2021.109118>)

**METHOD DETAILS, Single-cell RNA-seq analysis, page e4** - Analysis:
“To infer the lineage development of TEFF / EM to TCM, we used Monocle
package v2.5.4 (Qiu et al., 2017). The top 2000 significant DEGs were
chosen to order genes for trajectory reconstruction using the DDRTree
method followed by dimension reduction, cell trajectory inference and
pseudo-time measurements, which were computed via reversed graph
embedding.”

## Load data

Load Seurat objects into R

``` r
seulist <- readRDS(file.path(intermeddir, "seulist_drc.rds"))

trm <- seulist[["Trm"]]
disttrm <- seulist[["disTrm"]]
tem <- seulist[["Tem"]]
```

## 1. Merge Seurat objects

Merge the seurat objects together so that there is 1 object that
includes all samples.

Since we will be merging the cells, it is useful to add a prefix to the
cell names that identify which sample the cell came from. This
information can also be stored in the metadata.

``` r
Cells(trm)[1:10]
```

    ##  [1] "AAACCTGGTCTTGCGG" "AAACGGGAGGCTCAGA" "AAACGGGGTCAACATC" "AAACGGGTCTTGAGAC"
    ##  [5] "AAAGATGGTAGGGTAC" "AAAGATGTCCTACAGA" "AAAGCAACACACCGAC" "AAATGCCAGCTAGTTC"
    ##  [9] "AACCATGCACGGCCAT" "AACCATGCAGTTAACC"

``` r
trm.rn <- RenameCells(trm, add.cell.id = 'Trm')
Cells(trm.rn)[1:10]
```

    ##  [1] "Trm_AAACCTGGTCTTGCGG" "Trm_AAACGGGAGGCTCAGA" "Trm_AAACGGGGTCAACATC"
    ##  [4] "Trm_AAACGGGTCTTGAGAC" "Trm_AAAGATGGTAGGGTAC" "Trm_AAAGATGTCCTACAGA"
    ##  [7] "Trm_AAAGCAACACACCGAC" "Trm_AAATGCCAGCTAGTTC" "Trm_AACCATGCACGGCCAT"
    ## [10] "Trm_AACCATGCAGTTAACC"

``` r
# do this for the remaining samples that will be merged
tem.rn <- RenameCells(tem, add.cell.id='Tem')
disttrm.rn <- RenameCells(disttrm, add.cell.id='MucosalTrm')
```

Merge Trm and MucosalTrm cells

``` r
comb <- merge(trm.rn, disttrm.rn)

comb@meta.data %>%
  group_by(orig.ident) %>%
  summarize(ncells = length(orig.ident))
```

    ## # A tibble: 2 × 2
    ##   orig.ident ncells
    ##   <chr>       <int>
    ## 1 disTrm        863
    ## 2 Trm           459

Add the Tem cells

``` r
comb <- merge(comb, tem.rn)

comb@meta.data %>%
  group_by(orig.ident) %>%
  summarize(ncells = length(orig.ident))
```

    ## # A tibble: 3 × 2
    ##   orig.ident ncells
    ##   <chr>       <int>
    ## 1 disTrm        863
    ## 2 Tem          1457
    ## 3 Trm           459

Check the cell names

``` r
dim(comb@meta.data) #2779 cells in total
```

    ## [1] 2779    6

``` r
Cells(comb) %>%
  str_split("_") %>%
  map_dfr(~data.frame(prefix = .x[1], orig.cellname = .x[2])) %>%
  mutate(new.cellname = Cells(comb)) -> cellname.df

cellname.df %>%
  group_by(prefix) %>%
  summarize(ncells = length(orig.cellname))
```

    ## # A tibble: 3 × 2
    ##   prefix     ncells
    ##   <chr>       <int>
    ## 1 MucosalTrm    863
    ## 2 Tem          1457
    ## 3 Trm           459

## 2. Create a Monocle object

Monocle uses a different type of object, a “CellDataSet” object, to
store single-cell data. There are three data elements to a CellDataSet
object:

-   **cellData**: expression matrix with rows as genes and columns as
    cells
-   **phenoData**: data frame containing attributes of individual cells
-   **featureData**: data frame containing attributes of gene features

Extract the count matrix with rows as genes and columns as cells.

``` r
cd <- comb@assays$RNA@counts
dim(cd)
```

    ## [1] 13009  2779

``` r
rownames(cd)[1:5]
```

    ## [1] "Mrpl15"  "Lypla1"  "Tcea1"   "Atp6v1h" "Rb1cc1"

``` r
colnames(cd)[1:5]
```

    ## [1] "Trm_AAACCTGGTCTTGCGG" "Trm_AAACGGGAGGCTCAGA" "Trm_AAACGGGGTCAACATC"
    ## [4] "Trm_AAACGGGTCTTGAGAC" "Trm_AAAGATGGTAGGGTAC"

Extract the cell-level data.

``` r
pd <- comb@meta.data
dim(pd)
```

    ## [1] 2779    6

``` r
rownames(pd)[1:5]
```

    ## [1] "Trm_AAACCTGGTCTTGCGG" "Trm_AAACGGGAGGCTCAGA" "Trm_AAACGGGGTCAACATC"
    ## [4] "Trm_AAACGGGTCTTGAGAC" "Trm_AAAGATGGTAGGGTAC"

Extract the gene-level metadata.

``` r
rownames(comb)[1:5]
```

    ## [1] "Mrpl15"  "Lypla1"  "Tcea1"   "Atp6v1h" "Rb1cc1"

``` r
fd <- data.frame(gene_short_name = rownames(comb), 
                 row.names = rownames(comb))
dim(fd)
```

    ## [1] 13009     1

``` r
row.names(fd)[1:5]
```

    ## [1] "Mrpl15"  "Lypla1"  "Tcea1"   "Atp6v1h" "Rb1cc1"

Create the CellDataSet object.

``` r
pd.adf <- new("AnnotatedDataFrame", data = pd)
fd.adf <- new("AnnotatedDataFrame", data = fd)

cds <- newCellDataSet(cellData = cd, 
                      phenoData = pd.adf, 
                      featureData = fd.adf)
```

## 3. Run tranjectory analysis

Before running the trajectory analysis, the counts need to be
preprocessed in the following ways.

Estimate size factors. **Size factor corresponds with cell library size?
Is this similar to how size factors estimated for bulk DE in DESeq? Add
info?** Notice that the estimateSizeFactors() function calculates a size
factor for each cell and populates the column `Size_Factor` in the
phenoData dataframe.

``` r
cds.es <- estimateSizeFactors(cds)
head(pData(cds))
```

    ##                      orig.ident nCount_RNA nFeature_RNA percent.mt
    ## Trm_AAACCTGGTCTTGCGG        Trm       4480         1544   2.343750
    ## Trm_AAACGGGAGGCTCAGA        Trm       7628         1696   3.224961
    ## Trm_AAACGGGGTCAACATC        Trm       4830         1572   2.939959
    ## Trm_AAACGGGTCTTGAGAC        Trm       8058         2239   2.593696
    ## Trm_AAAGATGGTAGGGTAC        Trm       2378         1003   3.195963
    ## Trm_AAAGATGTCCTACAGA        Trm       2934         1027   2.453988
    ##                      RNA_snn_res.0.6 seurat_clusters Size_Factor
    ## Trm_AAACCTGGTCTTGCGG               1               1          NA
    ## Trm_AAACGGGAGGCTCAGA               0               0          NA
    ## Trm_AAACGGGGTCAACATC               0               0          NA
    ## Trm_AAACGGGTCTTGAGAC               1               1          NA
    ## Trm_AAAGATGGTAGGGTAC               0               0          NA
    ## Trm_AAAGATGTCCTACAGA               3               3          NA

``` r
head(pData(cds.es))
```

    ##                      orig.ident nCount_RNA nFeature_RNA percent.mt
    ## Trm_AAACCTGGTCTTGCGG        Trm       4480         1544   2.343750
    ## Trm_AAACGGGAGGCTCAGA        Trm       7628         1696   3.224961
    ## Trm_AAACGGGGTCAACATC        Trm       4830         1572   2.939959
    ## Trm_AAACGGGTCTTGAGAC        Trm       8058         2239   2.593696
    ## Trm_AAAGATGGTAGGGTAC        Trm       2378         1003   3.195963
    ## Trm_AAAGATGTCCTACAGA        Trm       2934         1027   2.453988
    ##                      RNA_snn_res.0.6 seurat_clusters Size_Factor
    ## Trm_AAACCTGGTCTTGCGG               1               1   1.2502888
    ## Trm_AAACGGGAGGCTCAGA               0               0   2.1288399
    ## Trm_AAACGGGGTCAACATC               0               0   1.3479676
    ## Trm_AAACGGGTCTTGAGAC               1               1   2.2488453
    ## Trm_AAAGATGGTAGGGTAC               0               0   0.6636577
    ## Trm_AAAGATGTCCTACAGA               3               3   0.8188275

Estimate dispersions **Add info about the purpose of this step ** Notice
that the estimated dispersion data can be extracted from the CellDataSet
object using the function dispersionTable().

``` r
cds.ds <- estimateDispersions(cds.es)

#head(dispersionTable(cds.es)) # if you haven't already estimated the dispersion, this will throw an error
head(dispersionTable(cds.ds))
```

    ##         gene_id mean_expression dispersion_fit dispersion_empirical
    ## 1        Mrpl15      0.22310305      1.0958020            0.0000000
    ## 2        Lypla1      0.14178132      1.5509077            0.0000000
    ## 3         Tcea1      0.31170279      0.8702657            0.2548530
    ## 4       Atp6v1h      0.13508195      1.6128301            0.1505939
    ## 5        Rb1cc1      0.14033596      1.5637671            0.3681927
    ## 6 4732440D04Rik      0.03058207      6.0908067            1.4588278

Filter out genes with low expression.

Notice that the detectGenes() function sets the expression detection
threshold, calculates the number of cells that express each gene above
that threshold, and adds that value to the feature dataframe.

``` r
head(fData(cds.ds))
```

    ##               gene_short_name
    ## Mrpl15                 Mrpl15
    ## Lypla1                 Lypla1
    ## Tcea1                   Tcea1
    ## Atp6v1h               Atp6v1h
    ## Rb1cc1                 Rb1cc1
    ## 4732440D04Rik   4732440D04Rik

``` r
cds.filt <- detectGenes(cds.ds, min_expr = 0.1) # set the detection threshold to 0.1
head(fData(cds.filt))
```

    ##               gene_short_name num_cells_expressed
    ## Mrpl15                 Mrpl15                 609
    ## Lypla1                 Lypla1                 410
    ## Tcea1                   Tcea1                 758
    ## Atp6v1h               Atp6v1h                 372
    ## Rb1cc1                 Rb1cc1                 373
    ## 4732440D04Rik   4732440D04Rik                  87

Identify genes that are expressed in at least 100 cells. These genes
will be retained.

``` r
dim(fData(cds.filt))[1] # 13009 features in total
```

    ## [1] 13009

``` r
fData(cds.filt) %>%
  filter(num_cells_expressed >= 100) %>%
  row.names() -> expressed_genes
expressed_genes[1:5]
```

    ## [1] "Mrpl15"  "Lypla1"  "Tcea1"   "Atp6v1h" "Rb1cc1"

``` r
length(expressed_genes) # 6543 features that are expressed in at least 100 cells
```

    ## [1] 6543

Select genes

``` r
cds.filt # assayData: 13009 features, 2779 samples
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 13009 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... num_genes_expressed (8 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Mrpl15 Lypla1 ... Adrb1 (13009 total)
    ##   fvarLabels: gene_short_name num_cells_expressed
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
cds.exprs <- cds.filt[expressed_genes,] 
cds.exprs # assayData: 6543 features, 2779 samples 
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 6543 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... num_genes_expressed (8 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Mrpl15 Lypla1 ... Wdr95 (6543 total)
    ##   fvarLabels: gene_short_name num_cells_expressed
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

Test if genes are differentially expressed based on sample type **Note
that I changed this to run using 2 cores because I think that is what
was requested by default for the workspace**

``` r
unique(pData(cds.exprs)$orig.ident)
```

    ## [1] "Trm"    "disTrm" "Tem"

``` r
#diff_test_res <- differentialGeneTest(cds.exprs, fullModelFormulaStr = "~orig.ident", cores = 2)
# since this can take some time, I saved the resulting object `diff_test_res` as an intermediate
#saveRDS(diff_test_res, file.path(intermeddir, "diff_test_res.RData"))

diff_test_res <- readRDS(file.path(intermeddir, "diff_test_res.RData"))
```

Examine DE results table. Tested genes are shown in the rows. Columns
include…

-   **status**: What does “OK” mean?
-   **family**: What does “negbinomial.size” mean?
-   **pval** : p value associated with gene
-   **qval**: q value, aka adjusted p value?
-   **gene_short_name** : Should match rownames
-   **num_cells_expressed** : Should match values in fData

``` r
dim(diff_test_res) # results table with rows as genes (n = 6543)
```

    ## [1] 6543    6

``` r
head(diff_test_res)
```

    ##         status           family        pval        qval gene_short_name
    ## Mrpl15      OK negbinomial.size 0.127561557 0.215779542          Mrpl15
    ## Lypla1      OK negbinomial.size 0.478860289 0.589054873          Lypla1
    ## Tcea1       OK negbinomial.size 0.135920966 0.227043880           Tcea1
    ## Atp6v1h     OK negbinomial.size 0.002470229 0.007849786         Atp6v1h
    ## Rb1cc1      OK negbinomial.size 0.331330990 0.444824392          Rb1cc1
    ## Pcmtd1      OK negbinomial.size 0.666407323 0.751646805          Pcmtd1
    ##         num_cells_expressed
    ## Mrpl15                  609
    ## Lypla1                  410
    ## Tcea1                   758
    ## Atp6v1h                 372
    ## Rb1cc1                  373
    ## Pcmtd1                  410

All have status == “OK” and family == “negbinomial.size”

``` r
diff_test_res %>%
  group_by(status, family) %>%
  summarize(n = length(gene_short_name))
```

    ## # A tibble: 1 × 3
    ## # Groups:   status [1]
    ##   status family               n
    ##   <chr>  <chr>            <int>
    ## 1 OK     negbinomial.size  6543

Pull out DEG based on qval \< 0.05. Use these genes to cluster/order
cells. Notice that after applying setOrderingFilter(), the column
“use_for_ordering” gets added to the featureData.

``` r
diff_test_res %>%
  filter(qval < 0.05) -> deg.df

dim(diff_test_res)[1] # 6543 genes tested
```

    ## [1] 6543

``` r
dim(deg.df)[1] # 2763 DEGs
```

    ## [1] 2763

``` r
ordering_genes <- row.names(subset(diff_test_res, qval < 0.05))
cds.o <- setOrderingFilter(cds.exprs, ordering_genes)

head(fData(cds.o)) # now has the column "use_for_ordering"
```

    ##         gene_short_name num_cells_expressed use_for_ordering
    ## Mrpl15           Mrpl15                 609            FALSE
    ## Lypla1           Lypla1                 410            FALSE
    ## Tcea1             Tcea1                 758            FALSE
    ## Atp6v1h         Atp6v1h                 372             TRUE
    ## Rb1cc1           Rb1cc1                 373            FALSE
    ## Pcmtd1           Pcmtd1                 410            FALSE

Reduce data dimensionality using Discriminative Dimensionality Reduction
with Trees (DDRTree).

From reduceDimension() help: “Prior to reducing the dimensionality of
the data, it usually helps to normalize it so that highly expressed or
highly variable genes don’t dominate the computation. The function
reduceDimension() automatically transforms the data in one of several
ways depending on the expressionFamily of the CellDataSet object. If the
expressionFamily is negbinomial or negbinomial.size, the data are
variance-stabilized.”

``` r
#cds.rd <- reduceDimension(cds.o, max_components = 2, method = 'DDRTree')
# since this can take some time, I saved the resulting object `cds.rd` as an intermediate
#saveRDS(cds.rd, file.path(intermeddir, "cds-rd.RData"))

cds.rd <- readRDS(file.path(intermeddir, "cds-rd.RData"))
cds.rd
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 6543 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... num_genes_expressed (8 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Mrpl15 Lypla1 ... Wdr95 (6543 total)
    ##   fvarLabels: gene_short_name num_cells_expressed use_for_ordering
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

Notice that data generated by reduceDimensions() are stored inside the
CellDataSet object

``` r
dim(cds.rd@reducedDimS)[2] # 2779 columns, coordinates for each cell
```

    ## [1] 2779

``` r
cds.rd@reducedDimS[1:2,1:5] # rows correspond to the dimensions, i.e. here we set max_components to 2
```

    ##      Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA Trm_AAACGGGGTCAACATC
    ## [1,]           5.71645658           -3.6653896            4.4251530
    ## [2,]          -0.03339826            0.4716709            0.3273042
    ##      Trm_AAACGGGTCTTGAGAC Trm_AAAGATGGTAGGGTAC
    ## [1,]             7.036108             3.199824
    ## [2,]            -1.546211             0.972328

Order cells along the trajectory.

From orderCells() help: “Learns a”trajectory” describing the biological
process the cells are going through, and calculates where each cell
falls within that trajectory. …This function takes as input a
CellDataSet and returns it with two new columns: Pseudotime and State,
which together encode where each cell maps to the trajectory in the
reduced-dimension space.

Notice that after applying orderCells(), the columns “Pseudotime” and
“State” get added to the phenoData.

``` r
# cds.oc <- orderCells(cds.rd)
# since this can take some time, I saved the resulting object `cds.oc` as an intermediate
#saveRDS(cds.oc, file.path(intermeddir, "cds-oc.RData"))

cds.oc <- readRDS(file.path(intermeddir, "cds-oc.RData"))
cds.oc
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 6543 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... State (10 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Mrpl15 Lypla1 ... Wdr95 (6543 total)
    ##   fvarLabels: gene_short_name num_cells_expressed use_for_ordering
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

``` r
head(pData(cds.oc)) # now has the columns "Pseudotime" and "State"
```

    ##                      orig.ident nCount_RNA nFeature_RNA percent.mt
    ## Trm_AAACCTGGTCTTGCGG        Trm       4480         1544   2.343750
    ## Trm_AAACGGGAGGCTCAGA        Trm       7628         1696   3.224961
    ## Trm_AAACGGGGTCAACATC        Trm       4830         1572   2.939959
    ## Trm_AAACGGGTCTTGAGAC        Trm       8058         2239   2.593696
    ## Trm_AAAGATGGTAGGGTAC        Trm       2378         1003   3.195963
    ## Trm_AAAGATGTCCTACAGA        Trm       2934         1027   2.453988
    ##                      RNA_snn_res.0.6 seurat_clusters Size_Factor
    ## Trm_AAACCTGGTCTTGCGG               1               1   1.2502888
    ## Trm_AAACGGGAGGCTCAGA               0               0   2.1288399
    ## Trm_AAACGGGGTCAACATC               0               0   1.3479676
    ## Trm_AAACGGGTCTTGAGAC               1               1   2.2488453
    ## Trm_AAAGATGGTAGGGTAC               0               0   0.6636577
    ## Trm_AAAGATGTCCTACAGA               3               3   0.8188275
    ##                      num_genes_expressed Pseudotime State
    ## Trm_AAACCTGGTCTTGCGG                1544   5.576941     3
    ## Trm_AAACGGGAGGCTCAGA                1696  15.386912     4
    ## Trm_AAACGGGGTCAACATC                1572   6.908313     3
    ## Trm_AAACGGGTCTTGAGAC                2239   3.603911     1
    ## Trm_AAAGATGGTAGGGTAC                1003   8.283603     3
    ## Trm_AAAGATGTCCTACAGA                1027  12.949319     4

Plot trajectory

``` r
plot_cell_trajectory(cds.oc, color_by = "orig.ident")
```

![](6_trajectory_files/figure-markdown_github/unnamed-chunk-22-1.png)

## 4. Test for DE over psuedotime

Pull out all “Il” genes. Subset the CellDataSet to the “Il” genes.

``` r
fData(cds.oc) %>%
  pull(gene_short_name) -> all.ordered.genes

marker_genes <- grep('^Il', all.ordered.genes, value = T)
marker_genes
```

    ##  [1] "Il18r1"  "Il18rap" "Ilkap"   "Il2ra"   "Il15ra"  "Il2rg"   "Il2"    
    ##  [8] "Il6ra"   "Ilf2"    "Il17ra"  "Il16"    "Ilk"     "Il4ra"   "Il21r"  
    ## [15] "Ilvbl"   "Il27ra"  "Ilf3"    "Il10ra"  "Il6st"   "Il7r"    "Il2rb"  
    ## [22] "Il10rb"

``` r
cds.mg <- cds.oc[marker_genes,]
cds.mg
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 22 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... State (10 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Il18r1 Il18rap ... Il10rb (22 total)
    ##   fvarLabels: gene_short_name num_cells_expressed use_for_ordering
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

Test whether the expression of each gene varies with pseudotime.

``` r
diff_test_res2 <- differentialGeneTest(cds.mg, fullModelFormulaStr = "~sm.ns(Pseudotime)")
diff_test_res2
```

    ##         status           family          pval          qval gene_short_name
    ## Il18r1      OK negbinomial.size 6.489983e-154 1.427796e-152          Il18r1
    ## Il18rap     OK negbinomial.size  6.257838e-13  2.753449e-12         Il18rap
    ## Ilkap       OK negbinomial.size  8.383425e-01  8.383425e-01           Ilkap
    ## Il2ra       OK negbinomial.size  4.315236e-22  3.164506e-21           Il2ra
    ## Il15ra      OK negbinomial.size  1.032385e-01  1.336027e-01          Il15ra
    ## Il2rg       OK negbinomial.size  8.718242e-06  1.743648e-05           Il2rg
    ## Il2         OK negbinomial.size  3.310799e-21  1.820940e-20             Il2
    ## Il6ra       OK negbinomial.size  3.876050e-08  1.065914e-07           Il6ra
    ## Ilf2        OK negbinomial.size  4.130785e-04  6.990559e-04            Ilf2
    ## Il17ra      OK negbinomial.size  9.314943e-10  2.927553e-09          Il17ra
    ## Il16        OK negbinomial.size  5.354651e-01  6.200122e-01            Il16
    ## Ilk         OK negbinomial.size  2.079073e-05  3.811634e-05             Ilk
    ## Il4ra       OK negbinomial.size  2.503606e-07  6.119926e-07           Il4ra
    ## Il21r       OK negbinomial.size  5.806064e-06  1.277334e-05           Il21r
    ## Ilvbl       OK negbinomial.size  6.454102e-01  7.099512e-01           Ilvbl
    ## Il27ra      OK negbinomial.size  9.117076e-03  1.337171e-02          Il27ra
    ## Ilf3        OK negbinomial.size  8.369843e-01  8.383425e-01            Ilf3
    ## Il10ra      OK negbinomial.size  1.126128e-02  1.548426e-02          Il10ra
    ## Il6st       OK negbinomial.size  2.914106e-26  3.205517e-25           Il6st
    ## Il7r        OK negbinomial.size  2.749496e-03  4.320637e-03            Il7r
    ## Il2rb       OK negbinomial.size  1.590325e-11  5.831191e-11           Il2rb
    ## Il10rb      OK negbinomial.size  3.932014e-01  4.805795e-01          Il10rb
    ##         num_cells_expressed use_for_ordering
    ## Il18r1                 1210             TRUE
    ## Il18rap                 216             TRUE
    ## Ilkap                   578            FALSE
    ## Il2ra                   161             TRUE
    ## Il15ra                  111            FALSE
    ## Il2rg                  1384            FALSE
    ## Il2                     461             TRUE
    ## Il6ra                   328             TRUE
    ## Ilf2                    524             TRUE
    ## Il17ra                  703             TRUE
    ## Il16                    817            FALSE
    ## Ilk                     593            FALSE
    ## Il4ra                   270             TRUE
    ## Il21r                   656             TRUE
    ## Ilvbl                   172            FALSE
    ## Il27ra                  634             TRUE
    ## Ilf3                    181             TRUE
    ## Il10ra                  140            FALSE
    ## Il6st                   259             TRUE
    ## Il7r                   1390            FALSE
    ## Il2rb                   793             TRUE
    ## Il10rb                  563            FALSE

Notice how the model formula is specified. The term “Pseudotime” refers
to the column in phenoData.

``` r
head(pData(cds.mg))
```

    ##                      orig.ident nCount_RNA nFeature_RNA percent.mt
    ## Trm_AAACCTGGTCTTGCGG        Trm       4480         1544   2.343750
    ## Trm_AAACGGGAGGCTCAGA        Trm       7628         1696   3.224961
    ## Trm_AAACGGGGTCAACATC        Trm       4830         1572   2.939959
    ## Trm_AAACGGGTCTTGAGAC        Trm       8058         2239   2.593696
    ## Trm_AAAGATGGTAGGGTAC        Trm       2378         1003   3.195963
    ## Trm_AAAGATGTCCTACAGA        Trm       2934         1027   2.453988
    ##                      RNA_snn_res.0.6 seurat_clusters Size_Factor
    ## Trm_AAACCTGGTCTTGCGG               1               1   1.2502888
    ## Trm_AAACGGGAGGCTCAGA               0               0   2.1288399
    ## Trm_AAACGGGGTCAACATC               0               0   1.3479676
    ## Trm_AAACGGGTCTTGAGAC               1               1   2.2488453
    ## Trm_AAAGATGGTAGGGTAC               0               0   0.6636577
    ## Trm_AAAGATGTCCTACAGA               3               3   0.8188275
    ##                      num_genes_expressed Pseudotime State
    ## Trm_AAACCTGGTCTTGCGG                1544   5.576941     3
    ## Trm_AAACGGGAGGCTCAGA                1696  15.386912     4
    ## Trm_AAACGGGGTCAACATC                1572   6.908313     3
    ## Trm_AAACGGGTCTTGAGAC                2239   3.603911     1
    ## Trm_AAAGATGGTAGGGTAC                1003   8.283603     3
    ## Trm_AAAGATGTCCTACAGA                1027  12.949319     4

Then, `Pseudotime` is surrounded by the function `sm.ns()`, which is
defined in the R package VGAM. It is a “smart” version of splines::ns()
that specifies the model as a spline regression.

``` r
?sm.ns()
?ns()
```

Select genes that are significant at an FDR \< 5%

``` r
diff_test_res2 %>%
  filter(qval < 0.05) -> sig_genes
sig_genes
```

    ##         status           family          pval          qval gene_short_name
    ## Il18r1      OK negbinomial.size 6.489983e-154 1.427796e-152          Il18r1
    ## Il18rap     OK negbinomial.size  6.257838e-13  2.753449e-12         Il18rap
    ## Il2ra       OK negbinomial.size  4.315236e-22  3.164506e-21           Il2ra
    ## Il2rg       OK negbinomial.size  8.718242e-06  1.743648e-05           Il2rg
    ## Il2         OK negbinomial.size  3.310799e-21  1.820940e-20             Il2
    ## Il6ra       OK negbinomial.size  3.876050e-08  1.065914e-07           Il6ra
    ## Ilf2        OK negbinomial.size  4.130785e-04  6.990559e-04            Ilf2
    ## Il17ra      OK negbinomial.size  9.314943e-10  2.927553e-09          Il17ra
    ## Ilk         OK negbinomial.size  2.079073e-05  3.811634e-05             Ilk
    ## Il4ra       OK negbinomial.size  2.503606e-07  6.119926e-07           Il4ra
    ## Il21r       OK negbinomial.size  5.806064e-06  1.277334e-05           Il21r
    ## Il27ra      OK negbinomial.size  9.117076e-03  1.337171e-02          Il27ra
    ## Il10ra      OK negbinomial.size  1.126128e-02  1.548426e-02          Il10ra
    ## Il6st       OK negbinomial.size  2.914106e-26  3.205517e-25           Il6st
    ## Il7r        OK negbinomial.size  2.749496e-03  4.320637e-03            Il7r
    ## Il2rb       OK negbinomial.size  1.590325e-11  5.831191e-11           Il2rb
    ##         num_cells_expressed use_for_ordering
    ## Il18r1                 1210             TRUE
    ## Il18rap                 216             TRUE
    ## Il2ra                   161             TRUE
    ## Il2rg                  1384            FALSE
    ## Il2                     461             TRUE
    ## Il6ra                   328             TRUE
    ## Ilf2                    524             TRUE
    ## Il17ra                  703             TRUE
    ## Ilk                     593            FALSE
    ## Il4ra                   270             TRUE
    ## Il21r                   656             TRUE
    ## Il27ra                  634             TRUE
    ## Il10ra                  140            FALSE
    ## Il6st                   259             TRUE
    ## Il7r                   1390            FALSE
    ## Il2rb                   793             TRUE

``` r
cds.mg1 <- cds.mg[sig_genes$gene_short_name,]
cds.mg1
```

    ## CellDataSet (storageMode: environment)
    ## assayData: 16 features, 2779 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: Trm_AAACCTGGTCTTGCGG Trm_AAACGGGAGGCTCAGA ...
    ##     Tem_TTTGTCACATCTCCCA (2779 total)
    ##   varLabels: orig.ident nCount_RNA ... State (10 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: Il18r1 Il18rap ... Il2rb (16 total)
    ##   fvarLabels: gene_short_name num_cells_expressed use_for_ordering
    ##   fvarMetadata: labelDescription
    ## experimentData: use 'experimentData(object)'
    ## Annotation:

Plot individual genes

``` r
plot_genes_in_pseudotime(cds.mg1, ncol=3)
```

![](6_trajectory_files/figure-markdown_github/unnamed-chunk-28-1.png)

## 5. Additional resources

Developers of Monocle, are releasing a new package version: Monocle3.
This still in beta phase, so use with caution.
<https://cole-trapnell-lab.github.io/monocle3/docs/updates/>
