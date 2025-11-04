
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

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(miRNARefine)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
