
# schnapps

<!-- badges: start -->
[![R build status](https://github.com/shahcompbio/schnapps/workflows/R-CMD-check/badge.svg)](https://github.com/shahcompbio/schnapps/actions)
[![Codecov test coverage](https://codecov.io/gh/shahcompbio/schnapps/branch/master/graph/badge.svg)](https://codecov.io/gh/shahcompbio/schnapps?branch=master)
[![Docker](https://img.shields.io/docker/cloud/build/marcjwilliams1/schnapps)](https://hub.docker.com/repository/docker/marcjwilliams1/schnapps)
<!-- badges: end -->

schnapps (Single Cell Haplotype copy Number Analysis through Phased Probabilistic States) is a tool to estimate allele and haplotype specific copy number states in single cells with low coverage (~0.01X). schnapps phases alleles based on losses and gains across all cells and then assigns allele specific states for each bin in each cell using a hidden markov model.  Documentation is available [here](https://shahcompbio.github.io/schnapps/).

You can read more about schnapps in our [preprint](https://www.biorxiv.org/content/10.1101/2021.06.04.447031v1).

## Installation

You can install schnapps with the following command 

``` r
devtools::install_github("shahcompbio/schnapps")
```

## Input data

`schnapps` was developed to work with Direct Library Preperation + (DLP+) data. A high throughput single cell whole genome sequencing workflow, described in [Laks et al.](https://www.sciencedirect.com/science/article/pii/S0092867419311766). As such it works using the output of the the pipeline developed to process this type of data, available [here](https://github.com/shahcompbio/single_cell_pipeline). Despite being developed with this type of data and pipeline in mind, it should work well with other single cell whole genome technologies. The required inputs are total copy number estimates in bins across the genome and haplotype block counts per cell (SNP counts may also work). See the test datasets provided with the package for example inputs. If you have a different type of technology and would like some advice or help running schnapps please open an issue.

Below is an example of how to use schnapps with DLP data. You will need the HMM copy results table (`CNbins`) with the following columns: `chr`, `start`,`end`, `cell_id`, `state`, `copy`, as well as cell specific haplotype counts as outputted by [scgenome](https://github.com/shahcompbio/scgenome) with the following command. This includes the following columns: `chr`, `start`,`end`, `cell_id`, `hap_label`, `allele_id`, `readcount`.

```py
allele_results = scgenome.loaders.allele.load_haplotype_allele_data(
    ticket_directory,
)
allele_data = allele_results['allele_counts']
```

## Example

First we need to do some data wrangling to convert the `allele_data` table from long to wide format.
``` r
library(schnapps)
allele_data <- format_haplotypes_dlp(allele_data, CNbins)
```

Then we can call haplotype specific state per cell:
```r
hscn <- callHaplotypeSpecificCN(CNbins, allele_data)
```

Or alternatively allele specific states:
```r
hscn <- callAlleleSpecificCN(CNbins, allele_data)
```

See the vignette for more information on the differences between these two outputs.

After performing this inference, to QC the results it is useful to plot a different representation of the BAF. Here we plot the BAF as a function of the inferred states. The black lines indicate where we should see the mean BAF based on the state.
``` r
plotBAFperstate(ascn, maxstate = 10)
```

After having ensured the results make sense, you can plot a heatmap of the states across all cells with the following.
```r
plotHeatmap(ascn, plotcol = "state_BAF")
```
This will cluster the cell using umap and hdbscan.

`schnapps` includes a number of other utilities for analysing single cell (haplotype-specific) copy number including the following:

* integration with scRNAseq using [seurat](https://satijalab.org/seurat/index.html)
* extensive plotting functions
* integration and plotting with structural variants
* clustering

Please see the vignettes for more information.





