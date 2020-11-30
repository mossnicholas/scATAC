Nucleosome Estimation
================

## Single-cell nucleosomes around TFBS

Below is an example clustering approach of nucleosome profiles of PBMCs
around the SMC3 TF motif. Here, it was found that six distinct clusters
appropriately characterize nucleosomal states around the SMC3. Different
clusters exhibit varying degrees of bimodal around the binding sites,
representing cells in different regulatory transition phases.

``` r
library(data.table)
library(RCurl)
source("https://raw.github.com/mossnicholas/scATAC/main/nucleosomes/scNuc_cluster.R")

PBMC_SMC = fread("https://raw.github.com/mossnicholas/scATAC/main/PBMC_mats/PBMC_SMC_cov_mat.txt")

ClusterNucs(coverage_matrix = PBMC_SMC, 
            clusters = 6, 
            center = FALSE)
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

We can generate heatmaps visualizing nucleosome profiles of single cells
that compose each clusters. The heatmaps are ordered by NFR
(nucleosomal-free region) score, a metric which assesses nucleosomal
bimodality of an individual cell.

Here is the heatmap corrosponding to cluster with the greatest
proportion of cells (33.8%). Notice the strong bimodality of cells in
this cluster

``` r
PBMC_cluster5 = fread("https://raw.github.com/mossnicholas/scATAC/main/PBMC_mats/PBMC_SMC_cluster5.txt")
HeatmapNucs(coverage_matrix = PBMC_cluster5,
            subset = T)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

Then the heatmap of the cluster with 20.5% of cells - there is
noticeable reduced definition of the +1 nucleosome. Prior to
transcription-factor binding, this nucleosome is often depleted as part
of an activation complex.

``` r
PBMC_cluster6 = fread("https://raw.github.com/mossnicholas/scATAC/main/PBMC_mats/PBMC_SMC_cluster6.txt")
HeatmapNucs(coverage_matrix = PBMC_cluster6,
            subset = F) 
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->
