# isoviz <img src="inst/figures/isoviz.png" align="right" height="139"/>

<!-- badges: start -->

[![GitHub issues](https://img.shields.io/github/issues/daklab/isoviz)](https://img.shields.io/github/issues/daklab/isoviz/issues) [![Travis build status](https://travis-ci.com/karini925/isoviz.svg?branch=master)](https://travis-ci.com/karini925/isoviz)

<!-- badges: end -->

`isoviz` is a package that allows for simplified transcript isoforms visualizations.

The goal of isoviz is to simplify working with trancsript isoforms. Given a gene name or Ensembl ID, you can visualize all transcript structures. Additionally, Isoviz integrates with junction coordinates and counts obtained with Regtools. You can choose to visualize only the transcript isoforms that are detected in your sample and see how they group into leafuctter intron clusters.

## Installation

Note, the package has thus far been tested under R version 4.3.0.

You can install the development version of isoviz from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("daklab/isoviz")
```

Alternatively, you can do the following:

``` r
# 1. git clone https://github.com/daklab/isoviz.git
# 2. set the directory in which isoviz is cloned into as your working directory, then do the following:
library(devtools)
devtools::build()
devtools::install()
library(isoviz)
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

``` r
exon_coords <- all_coordinates[[1]]
```

#### Let's look at the intron coordinates

``` r
intron_coords <- all_coordinates[[2]]
```

#### Let's make a basic plot of the isoforms for the gene RBFOX2

``` r
# get exons and introns 
rbfox2_exons <- filter(exon_coords, gene_name == "RBFOX2")
rbfox2_introns <- filter(intron_coords, gene_name == "RBFOX2")
```

#### Intron re-scaling

Let's rescale the intron coordinates but keep the relative exon alignments this step will return an updated dataframe of exon coordinates that we can use with our plotting function.

``` r
rescaled_coords = isoviz_rescale_introns(rbfox2_introns, rbfox2_exons, width_rescale=10) 
```

#### Intron clustering using junction count data from BAM files

We will first need to load our leafcutter junctions that we obtained by running Regtools extract junctions on our BAM files. You can use one of our preloaded cell types or first run this on your own BAM file and then use as input for this function. We will look at junctions in hESC data. We will first need to run 'minicutter' to cluster the junction coordinates and obtain intron cluster events.

``` r

# load junctions for cell type of interest or input your own
junctions <- system.file("data", "21c_hESC_RB2_empty_v41_basic.junc", package="isoviz")

# run minicutter to get clusters 
intron_clusts <- isoviz_minicutter(juncs_file = junctions)
```

Let's take a look at how junctions were grouped into intron clusters. Note, some of these clusters will only contain one junction (singleton). You will have the option to change this when running isoviz_minicutter.

``` r
# For now we will use the expanded intron dataset with additional annotations by Megan 
# We will need an additional function to add these annotations 
intron_annotations <- system.file("data", "gencode_intron_all_data.rda", package="isoviz")
load(intron_annotations)
```

#### Plotting isoforms and junction and cluster information

Let's continue working on RBFOX2. Now we have to map the observed junctions with their corresponding exons and transcript isoforms. To do this, we will use the `isoviz_map_junctions` function. Make sure to check strand of gene!

``` r
mapped_junctions = isoviz_map_junctions(cell_type = "hESC", gene_introns, intron_clusts, gencode_intron_all_data)
```

Now we are ready to make a plot!

``` r
isoviz_plot_juncs_to_iso(mapped_junctions, gene_exons, gene_introns,
                                    cell_type = "hESC",
                                    junc_usage = 5, #min junc usage to be included 
                                    include_all_juncs = TRUE, intron_scale = "no")
```

<img src="inst/figures/firstreadmeimage.png"/>

``` r
isoviz_plot_juncs_to_iso(mapped_junctions, gene_exons, gene_introns,
                                    cell_type = "hESC",
                                    junc_usage = 5, #min junc usage to be included 
                                    intron_scale = "yes", intron_scale_width = 10,
                                    include_all_juncs = TRUE)
```

<img src="inst/figures/secondreadmeimage.png"/>

``` r
isoviz_plot_juncs_to_iso(mapped_junctions, gene_exons, gene_introns,
                                    cell_type = "hESC",
                                    junc_usage = 5, #min junc usage to be included 
                                    intron_scale = "yes", intron_scale_width = 10,
                                    include_all_juncs = FALSE,
                                    include_specific_isoforms = c("RBFOX2-209", "RBFOX2-220", "RBFOX2-208", "RBFOX2-205"),
                                    include_specific_junctions = c("junc178147", "junc178149", "junc178135", "junc178136", "junc178145", "junc178146"))
```

<img src="inst/figures/thirdreadmeimage.png"/>

Finally, if you are interested in obtaining a list of sgRNAs that you could use to target specific junctions with Cas13, you can do the following:

``` r
guide_table = isoviz_get_guide_predictions(gene = "ENSG00000100320", 
                                           leafcutter_input=intron_clusts)
guide_table
```

## Citation

``` r
citation("isoviz")
#> 
#> To cite package ‘isoviz’ in publications use:
#>
#>  Isaev K, Schertzer M (2023). _isoviz: Visualize transcript isoforms alongside exon-exon junction counts_. R package version 0.1.0.
#>
#> A BibTeX entry for LaTeX users is
#>
#>  @Manual{,
#>    title = {isoviz: Visualize transcript isoforms alongside exon-exon junction counts},
#>    author = {Karin Isaev and Megan Schertzer},
#>    year = {2023},
#>    note = {R package version 0.1.0},
#>  }
```

## Questions or suggestions?
