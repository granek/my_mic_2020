# README bioinformatics

## 1_cellranger-demo.Rmd

**Goal of this workshop:** Learn how to preprocess scRNAseq data using cellranger

**What's covered in this workshop:**

- Transform BCL files into fastq files
- Map reads and count genes

## 2_qc.Rmd

**Goal of this workshop:** Learn how to QC scRNAseq data

**What's covered in this workshop:**

- Convert cellranger output to a Seurat object
- Calculate mitrochondiral gene percentage
- Learn data quality attributes
- Exclude data from low-quality cells

## 3_load_transformCounts.Rmd

**Goal of this workshop:** Learn how to load cellranger output into R and prepare QC'd data for downstream analyses

**What's covered in this workshop:**

- Load data into R as a Seurat object
- Verify that the data have been QC'd
- Transform gene counts

## 4_clustering_annotation.Rmd

**Goal of this workshop:** Illustrate the clustering and cell-type annotation procedures. 

**What's covered in this workshop:**

- Highly variable gene identification
- Dimension reduction
- Clustering and visualization
- Cell-type specific marker gene detection
- How to generate the figures for the paper?

## 5_de.Rmd

**Goal of this workshop:** Learn how to test differential expression between cells

**What's covered in this workshop:**

- Identify DEGs between clusters
- Re-name clusters
- Verify manuscript DE results

## 6_trajectory.Rmd

**Goal of this workshop:** Learn to analyze cell trajectories with monocle

**What's covered in this workshop:**

- Merge Seurat objects
- Create a Monocle object
- Run trajectory analysis
- Test for DE over pseudotime

---

# Data sources

The scripts `cellranger-demo.Rmd` and `qc.Rmd` use scRNA-seq dataset of 8K human peripheral blood mononuclear cells (PBMCs) freely available from 10X Genomics

- https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz\
- https://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
- http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_fastqs.tar

All downstream scripts use data from Christian et al. 2021 Cell Reports (https://doi.org/10.1016/j.celrep.2021.109118) 

