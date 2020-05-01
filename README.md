
# schnapps

schnapps (Single Cell Haplotype copy Number Analysis through Phased Probabilistic States) is a tool to estimate allele and haplotype specific copy number states in single cells with low coverage (~0.01X). schnapps phases alleles based on losses and gains across all cells and then assigns allele specific states for each bin in each cell using a hidden markov model. 

## Installation

You can install schnapps with the following command 

``` r
devtools::install_github("shahcompbio/schnapps")
```

## Input data

Below is an example of how to use schnapps with DLP data. You will need the HMM copy results table (`CNbins`) with the following columns: `chr`, `start`,`end`, `cell_id`, `state`, `copy`, as well as cell specific haplotype counts as outputted by [scgenome](https://github.com/shahcompbio/scgenome) with the following command.

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

Then we can call the allele specific states:
```r
ascn <- callAlleleSpecificCN(CNbins, allele_data)
```

Or alterantively the haplotype specific states:
```r
hscn <- callHaplotypeSpecificCN(CNbins, allele_data)
```

See the vignette for more information on the differences between these two outputs.

After performing this inference, to QC the results it is useful to plot a different representation of the BAF. Here we plot the BAF as a function of the inferred states. The black lines indicate where we should see the mean BAF based on the state.
``` r
plotBAFperstate(ascn, maxstate = 10)
```

After having ensured the results make sense, you can plot a heatmap of the states across all cells with the following.
```r
plotHeatmap(alleleCN, plotcol = "state_BAF")
```
This will cluster the cell using umap and hdbscan.




