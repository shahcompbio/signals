
# schnapps

<!-- badges: start -->
[![R build status](https://github.com/shahcompbio/schnapps/workflows/R-CMD-check/badge.svg)](https://github.com/shahcompbio/schnapps/actions)
[![Codecov test coverage](https://codecov.io/gh/shahcompbio/schnapps/branch/master/graph/badge.svg)](https://codecov.io/gh/shahcompbio/schnapps?branch=master)
[![Docker](https://img.shields.io/docker/cloud/build/marcjwilliams1/schnapps)](https://hub.docker.com/repository/docker/marcjwilliams1/schnapps)
<!-- badges: end -->

schnapps (Single Cell Haplotype copy Number Analysis through Phased Probabilistic States) is a tool to estimate allele and haplotype specific copy number states in single cells with low coverage (0.01-0.1X). schnapps phases alleles based on losses and gains across all cells and then assigns allele specific states for each bin in each cell using a hidden markov model.  Documentation is available [here](https://shahcompbio.github.io/schnapps/).

You can read more about schnapps in our [preprint](https://www.biorxiv.org/content/10.1101/2021.06.04.447031v1).

## Installation

You can install schnapps with the following command 

``` r
devtools::install_github("shahcompbio/schnapps")
```

## Input data

`schnapps` was developed to work with Direct Library Preperation + (DLP+) data. A high throughput single cell whole genome sequencing workflow, described in [Laks et al.](https://www.sciencedirect.com/science/article/pii/S0092867419311766). As such it works using the output of the the pipeline developed to process this type of data, available [here](https://github.com/shahcompbio/single_cell_pipeline). Despite being developed with this type of data and pipeline in mind, it should work well with other single cell whole genome technologies. We have run it succesfully with 10X CNV data for example. The required inputs are total copy number estimates in bins across the genome and haplotype block counts per cell (SNP counts may also work). See the test datasets provided with the package for example inputs. If you have a different type of technology and would like some advice or help running schnapps please open an issue. We describe in more detail the necessary input below.

### DLP+ data

You will need the HMM copy results table (`CNbins`) with the following columns: `chr`, `start`,`end`, `cell_id`, `state`, `copy`. `state` is the inferred total copy number state. `copy` values are GC-correceted, ploidy corrected normalized read counts (this is what HMMcopy uses to infer the total copy number states). You will also need the cell specific haplotype counts (`allele_data`) as outputted by the `inferhaps` and `couthaps` sub commands in the pipeline. This includes the following columns: `chr`, `start`,`end`, `cell_id`, `hap_label`, `allele_id`, `readcount`. `allele_id` gives the haplotype phasing, often this is designated as `1|0` or `0|1`, you'll notice in the test data each haplotype has 2 rows, one for each phase. `hap_label` is a label given to the haplotype block, these labels are consistent across cells, and in combination with `chr` gives each block a unique ID. These blocks are computed by finding SNP's that can be confidently phased together. See the section of "phasing uncertainty" in SHAPEIT [here](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#uncertainty). This notion of phasing uncertainty is also used in [remixt](https://github.com/amcpherson/remixt). Note that this notion of haplotype block is different to other tools such as CHISEL where it is assumed that the the phasing remains consistent over seem specified length such as 50kb. Finally the `readcount` columns gives the sum of the number of reads for each allele in each block

### Other technologies

Other technologies and software should also be compatible with schnapps. For example, we have used 10X data successfully. If you have single cell bam files or fastq files see the detailed documentation for running our single cell pipeline [here](https://github.com/shahcompbio/single_cell_pipeline/blob/master/docs/source/install.md). Alternatively, we provide a lightweight snakemake pipeline with the key steps [here](https://github.com/marcjwilliams1/hscn_pipeline). Also included there are some scripts to demultiplex 10X CNV bams.

If you total copy number calls from other software such as 10X cellranger or [scope](https://github.com/rujinwang/SCOPE), these may also work but not something we have tried. Feel free to open an issue if you need some advice.

### Data conventions

Some of the plotting tools assume that the `cell_id`'s conform to the following naming conventions

```
{sample}-{library}-{R*}-{C*}
```

Here R and C refer the rows and columns on the chip. If you're using another technology and your cells are named differently, we would reccomend renaming your cell's for easy compatability. For example, if you have 10X data where cell_id's are barcodes that look like `CCGTACTTCACGGTAT-1` something like this would work

```{r}
new_cell_id <- paste("mysample", "mylibrary", "CCGTACTTCACGGTAT-1", sep = "-")
```

It is imortant to have 4 string's seperated by "-", but the unique cell identifier should be one of the last 2 strings for the heatmap labelling to format nicely.

## Example

Here we'll show an example of running schnapps with a small dataset of 250 cells from the DLP platform. First we need to do some data wrangling to convert the `allele_data` table from long to wide format.

``` r
library(schnapps)
haplotypes<- format_haplotypes_dlp(haplotypes, CNbins)
```

Then we can call haplotype specific state per cell:
```r
hscn <- callHaplotypeSpecificCN(CNbins, haplotypes)
```

Or alternatively allele specific states:
```r
hscn <- callAlleleSpecificCN(CNbins, haplotypes)
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
* utilities such as consensus copy number in cell clusters and arm level changes etc

Please see the vignettes for more information.





