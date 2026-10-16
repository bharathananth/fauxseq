
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fauxseq

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

A Synthetic Timeseries RNA-seq Data Generator

## How to install fauxseq

You can install the development version of fauxseq from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("bharathananth/fauxseq")
```

## What fauxseq does

This package can generate synthetic timeseries RNA-seq data with a
fraction of rhythmic genes and empirical relationship between mean
expression and variability across replicates of genes.

`fauxseq` improves upon previous tools as follows: - generates
negative-binomial count data with empirically-estimated mean-dispersion
properties. - generates data with either differential expression or
differential rhythmicity or both.

## Get started with fauxseq

``` r
library(fauxseq)
sim_data <- generate_rhythmic_rnaseq()
head(sim_data$counts) # count data
#>     samples
#> id   qvsfc jsrni hcasm drgmf nquhy eaxzn
#>   g1   937  1171   959  1561   767   710
#>   g2    23    19    20    14    18    23
#>   g3  1707  1641  1332  1358  1490  1663
#>   g4  6391 92524  4395  7079 63037 19629
#>   g5  4986  8981  6393  8518  9145  9550
#>   g6  1666  1645  1469  1246  1728  1901
head(sim_data$params) # parameters of the rhythmic genes
#>   id         A      phi           mu  dispersion
#> 1 g1 0.0000000 0.000000   430.333333 0.171881593
#> 2 g2 0.0000000 0.000000     8.555556 0.011100650
#> 3 g3 0.0000000 0.000000   717.000000 0.013414883
#> 4 g4 0.2441755 4.601917 17156.000000 1.570429065
#> 5 g5 0.0000000 0.000000  3705.555556 0.044844613
#> 6 g6 0.0000000 0.000000   742.000000 0.009458109
head(sim_data$exp_design) # time labels of each sample (columns) of count data
#>       time
#> qvsfc    0
#> jsrni    4
#> hcasm    8
#> drgmf   12
#> nquhy   16
#> eaxzn   20
```

## How to cite fauxseq

Please cite this software using CITATION.CFF or the “Cite this
repository” link in the right sidebar.
