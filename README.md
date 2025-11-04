
<!-- README.md is generated from README.Rmd. Please edit that file -->

# miRNARefine

<!-- badges: start -->
<!-- badges: end -->

Useful Tool to Refine and Improve the Quality and Interpretability of
miRNA Expression Data.

## Description

`miRNARefine` is an R package which aims to identify and correct
potential data irregularities such as outliers, batch effects, and
instability in miRNA expression across samples, and to provide
diagnostic visualizations that help researchers assess data reliability
before downstream analyses.

`miRNARefine` was developed in the following environment:  
R version: 4.5.1 (2025-06-13) - “Great Square Root”  
Platform: Windows 11 x64(x86_64, mingw32)

## Installation

You can install the development version of miRNARefine like so:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("feiyangsun246/miRNARefine", build_vignettes = TRUE)
library("TestingPackage")
```

To run the Shiny app:

``` r
# Under construction
```

## Overview

To list all functions available in the package:

``` r
ls("package:miRNARefine")
```

`miRNARefine` contains 7 functions:

1.  **adaptiveFiltering** Automatically filters miRNAs based on
    expression, variance, or missing values.
2.  **detectOutliersPCA** Identifies anomalous samples using principal
    component analysis.
3.  **missingValueHandling** Detects and imputes missing values in the
    miRNA dataset. Supports adaptive methods like mean, median, or KNN
    imputation.
4.  **compareNormalization** Applies multiple normalization methods
    (e.g., log2, z-score, quantile) and compares their effects.
5.  **detectBatch** Detects and optionally corrects batch effects in the
    dataset.
6.  **miRNAStability** Calculates feature-level stability metrics such
    as coefficient of variation (CV) or median absolute deviation (MAD).
7.  **plotStabilityDistribution** Visualizes the distribution of miRNA
    stability metrics across all features. Highlights highly stable or
    highly variable miRNAs, helping users quickly assess feature
    quality.

<br>

To view datasets available in the package:

``` r
data(package = "miRNARefine") 
```

<br>

To view the tutorial of the package:

``` r
browseVignettes("miRNARefine")
```

The package also contains two datasets, called miRNASeq1 and miRNASeq2.
Refer to package vignettes for more details. An overview of the package
workflow is illustrated below.

![](inst/extdata/workflow.png)

## Contributions

The author of the package is Feiyang Sun. Feiyang designed the overall
package structure, implemented all functions, wrote the documentation
and vignettes, and performed testing and debugging. All analyses,
visualizations, and examples included in the package were created by
Feiyang, who is also responsible for maintaining and updating the
package. `preprocessCore` is used to Provide quantile normalization
methods for miRNA expression data. `stats` is used for feature stability
and outlier detection. `impute` provides functions for KNN-based
imputation of missing values. `sva` helps implement batch effect
detection and correction methods. `ggplot2` is used for creating
stability distribution plots.

`ChatGPT-5` from OpenAI was used to polish and format documentation
structure, debug code, and look up function information.

## References
