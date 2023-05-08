# isoviz <img src="inst/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/daklab/isoviz)](https://img.shields.io/github/issues/daklab/isoviz/issues)
<!-- badges: end -->

`isoviz` is a package that allows for simplified transcript isoforms visualizations.

The goal of isoviz is to simplify working with trancsript isoforms. Given a gene name or Ensembl ID, you can visualize all transcript structures. Additionally, Isoviz integrates with junction coordinates and counts obtained with Regtools. You can choose to visualize only the transcript isoforms that are detected in your sample and see how they group into leafuctter intron clusters. 

## Installation

You can install the development version of isoviz from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("daklab/isoviz")
```

## Introduction 

Our package aims to help quantify and visualize the transcript isoforms that are present in your samples.

## Example

We start off by loading all our exon and intron coordinates. By default, this will use 'gencode.v41.basic annotation' data. You can also input your own psl file (gtf file as well in the future).  

``` r
library(isoviz)

## basic example code
file_path <- system.file("data", "gencode.v41.basic.annotation.psl", package="isoviz")
gene_trans <- system.file("data", "gencode_v41_gene_transcript_convert.txt", package="isoviz")

all_coordinates <- isoviz_coords(file_path, gene_trans, input_type="psl") #use default genome .psl file  
```

#### Let's look at the exon coordinates 

```{r}
exon_coords <- all_coordinates[[1]]
print(head(exon_coords))
``` 

#### Let's look at the intron coordinates 

```{r}
intron_coords <- all_coordinates[[2]]
print(head(intron_coords))
``` 


## Questions or suggestions? 

