# Tutorial of single-cell RNA-seq data analysis in R
#### Compiled by Zhisong He, Barbara Treutlein
#### Updated on 2024-04-18
### Table of Content
  * [Introduction](#introduction)
  * [Preparation](#preparation)
  * [Now let's start Part 1](#now-lets-start-part-1)
    * [Step 0. Import Seurat package](#step-0-import-seurat-package)
    * [Step 1. Create a Seurat object](#step-1-create-a-seurat-object)
    * [Step 2. Quality control](#step-2-quality-control)
    * [Step 3. Normalization](#step-3-normalization)
    * [Step 4. Feature selection for following heterogeneity analysis](#step-4-feature-selection-for-following-heterogeneity-analysis)
    * [Step 5. Data scaling](#step-5-data-scaling)
    * [(Optional and advanced) Alternative step 3-5: to use SCTransform](#optional-and-advanced-alternative-step-3-5-to-use-sctransform)
    * [Step 6. Linear dimension reduction using principal component analysis (PCA)](#step-6-linear-dimension-reduction-using-principal-component-analysis-pca)
    * [Step 7. Non-linear dimension reduction for visualization](#step-7-non-linear-dimension-reduction-for-visualization)
    * [Step 8. Cluster the cells](#step-8-cluster-the-cells)
    * [Step 9. Annotate cell clusters](#step-9-annotate-cell-clusters)
    * [Step 10. Pseudotemporal cell ordering](#step-10-pseudotemporal-cell-ordering)
    * [Step 11. Save the result](#step-11-save-the-result)
    * [What else?](#what-else)
  * [Now starts Part 2: when you need to jointly analyze multiple scRNA-seq data sets](#now-starts-part-2-when-you-need-to-jointly-analyze-multiple-scrna-seq-data-sets)
    * [Step 0. Load data](#step-0-load-data)
    * [Step 1. Merge the two data sets](#step-1-merge-the-two-data-sets)
    * [Step 2-1. Data integration using Seurat](#step-2-1-data-integration-using-seurat)
    * [Step 2-2. Data integration using Harmony](#step-2-2-data-integration-using-harmony)
    * [Step 2-3. Data integration using LIGER](#step-2-3-data-integration-using-liger)
    * [Step 2-4. Data integration using MNN](#step-2-4-data-integration-using-mnn)
    * [Step 2-5. Data integration using RSS to BrainSpan](#step-2-5-data-integration-using-rss-to-brainspan)
    * [Step 2-6. Data integration using CSS](#step-2-6-data-integration-using-css)
    * [Step 3. How shall we compare different data integration methods](#step-3-how-shall-we-compare-different-data-integration-methods)
  * [Now starts Part 3: when you have an annotated reference data set and want it to facilitate the analysis of a new data](#now-starts-part-3-when-you-have-an-annotated-reference-data-set-and-want-it-to-facilitate-the-analysis-of-a-new-data)
    * [Step 0. Load data](#step-0-load-data)
    * [Method 1-1. Transcriptome similarity on cell cluster level](#method-1-1-transcriptome-similarity-on-cell-cluster-level)
    * [Method 1-2. Transcriptome similarity on cell level](#method-1-2-transcriptome-similarity-on-cell-level)
    * [Method 2. Seurat-based label transfer](#method-2-seurat-based-label-transfer)
    * [Other methods, and more to say](#other-methods-and-more-to-say)
  * [Now starts Part 4: more optional advanced analysis for scRNA-seq data](#now-starts-part-4-more-optional-advanced-analysis-for-scrna-seq-data)
    * [Part 4-1. Cluster connectivity analysis with PAGA](#part-4-1-cluster-connectivity-analysis-with-paga)
    * [Part 4-2. Pseudotime reconstruction without subseting into an unbranched trajectory](#part-4-2-pseudotime-reconstruction-without-subseting-into-an-unbranched-trajectory)
    * [Part 4-3. RNA velocity analysis](#part-4-3-rna-velocity-analysis)
    * [Part 4-4. Trajectory analysis with CellRank](#part-4-4-trajectory-analysis-with-cellrank)
    * [Part 4-5. Cell communication analysis](#part-4-5-cell-communication-analysis)


## Introduction
After getting the scRNA-seq data of your samples, you will want to analyze it properly.

Multiple toolkits and analytic frameworks have been developed to facilitate scRNA-seq data analysis. These options include but are not limit to [Seurat](https://satijalab.org/seurat/), developed by Rahul Satija's Lab in R, and [scanpy](https://icb-scanpy.readthedocs-hosted.com/en/stable/), developed by Fabian Theis's Lab in Python. Both toolkits provide functions and rich parameter sets that serve most of the routine analysis that one usually does on scRNA-seq data. However, one should be aware that these analytic frameworks do not cover all the interesting analyses that one can do when analyzing data. It is also important to get to know other tools for scRNA-seq data analysis.

Since this is a tutorial for beginners, we will mostly introduce how to use Seurat to analyze your scRNA-seq data in R. At the end, we will also mention some other additional tools (e.g. presto, destiny, Harmony, simspec, etc.), which provide additional functionalities that you may miss if you only use Seurat. In the most recent update, we also provide the briefly example of some commonly used advanced analysis, such as RNA velocity.

## Preparation
This tutorial assumes that the sequencing data preprocessing steps, including base calling, mapping and read counting, have been done. 10x Genomics has its own analysis pipeline [Cell Ranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) for data generated with the 10x Genomics Chromium Single Cell Gene Expression Solution. At the end of the Cell Ranger pipeline, a count matrix is generated. If your scRNA-seq data is generated using another technology (e.g. well-based experiments using Smart-Seq2 and others), the Cell Ranger pipeline is likely unapplicable, and you will have to find another solution to generate the count matrix.

As part of this tutorial, we are providing two data sets (DS1 and DS2), both generated using 10x Genomics and preprocessed using Cell Ranger. They are both public scRNA-seq data of human cerebral organoids and are part of the data presented in this [paper](https://www.nature.com/articles/s41586-019-1654-9). The first part of this tutorial, which includes most of the general analysis pipeline, is based on DS1, while the second part, which focuses on data integration and batch effect correction, is based on both data sets.

As a test for yourself, please try to apply what is learned in the first part to DS2 and only then continue with part 2 of the vignette. This will also give you an idea which types of cells are in DS2 and how comparable it is to DS1, before doing any data integration of both data sets.

## Now let's start Part 1
### Step 0. Import Seurat package
First of all, please make sure that Seurat is installed in your R.
```R
library(Seurat)
```
This imports your installed Seurat package into your current R session. No error should be seen but some verbose information is likely. If it warns you that the package is unavailable, please install Seurat first
```R
install.packages("Seurat")
library(Seurat)
```
### Step 1. Create a Seurat object
Seurat implements a new data type which is named 'Seurat'. It allows Seurat to store all the steps and results along the whole analysis. Therefore, the first step is to read in the data and create a Seurat object. Seurat has an easy solution for data generated using the 10x Genomics platform.
```R
counts <- Read10X(data.dir = "data/DS1/")
seurat <- CreateSeuratObject(counts, project="DS1")
```
What the ```Read10X``` function does is to read in the matrix and rename its row names and col names by gene symbols and cell barcodes, respectively. Alternatively, one can do this manually, which is probably what one would do when the data is not generated using 10x.
```R
library(Matrix)
counts <- readMM("data/DS1/matrix.mtx.gz")
barcodes <- read.table("data/DS1/barcodes.tsv.gz", stringsAsFactors=F)[,1]
features <- read.csv("data/DS1/features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
rownames(counts) <- make.unique(features[,2])
colnames(counts) <- barcodes

seurat <- CreateSeuratObject(counts, project="DS1")
```
If you look at the [Seurat tutorial](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html), you would notice that some extra options are added to the ```CreateSeuratObj``` function, such as ```min.cells``` and ```min.features```. When these two parameters are set, an initial filtering is applied to the data, removing right from the beginning all genes with reads detected in too few cells, as well as cells with too few genes detected. This is fine, but I personally recommend to keep all genes (i.e. default or ```min.cells = 0```)

### Step 2. Quality control
After creating the Seurat object, the next step is to do quality control on the data. The most common quality control is to filter out
1. Cells with too few genes detected. They usually represent cells which are not sequenced deep enough for reliable characterization.
2. Cells with too many genes detected. They may represent doublets or multiplets (i.e. two or more cells in the same droplet, therefore sharing the same cell barcode).
3. Cells with high mitochondrial transcript percentage. As most of the scRNA-seq experiments use oligo-T to capture mRNAs, mitochondrial transcripts should be relatively under-representative due to their lack of poly-A tails, but it is unavoidable that some mitochondrial transcripts are captured. Meanwhile, there is also some evidence that stable poly-A tails exist in some mitochondrial transcripts but serve as a marker for degradation (e.g. this [paper](https://mcb.asm.org/content/25/15/6427.long)). Together, cells with high mitochondrial transcript percentage likely represent cells under stress (e.g. hypoxia) which produce a lot of mitochondria, or which produce an abnormally high amount of truncated mitochondrial transcripts.

While numbers of detected genes are summarized by Seurat automatically when creating the Seurat object (with nFeature_RNA being the number of detected genes/features; nCount_RNA being the number of detected transcripts), one needs to calculate mitochondial transcript percentages manually. Still, Seurat provides an easy solution
```R
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
```

Please note that there is no one-size-fits-all filtering criteria, as the normal ranges of these metrics can vary dramatically from one experiment to another, depending on sample origin as well as reagents and sequencing depths. One suggestion here is to **ONLY FILTER OUT OUTLIER CELLS**, i.e. the **minority** of cells with certain QC metrics clearly above or below the majority of cells. To do that, one needs to first know how these values are distributed in the data. One can look at the distribution by creating a violin plot for each of the metrics.
```R
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
<img src="images/vlnplot_QC.png" align="centre" /><br/><br/>
Or if you don't like the dots (individual cells)
```R
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
```
