# signals 0.14.1

* Fix rephasebins stability to ensure idempotent results

# signals 0.14.0

* Add gene annotation feature to plotHeatmap with customizable positions and styling
* Add chromosome ideogram (cytoband) visualization at the bottom of heatmaps
* Add plotideogram and plotallbins parameters to visualize centromeric regions
* Improve frequency annotation handling for BAF and state plots

# signals 0.13.1

* Fix plotHeatmap annotation handling for data.table/tibble inputs

# signals 0.13.0

* Add SV visualization with lines and arcs style showing orientation
* Add automatic position and strand correction for SV data with reversed coordinates
* Add SV legend support for orientation color coding
* Fix SV read count axis scaling issues

# signals 0.12.0

* option to input cells to use for phasing for each chromosome
* plotHeatmap with arbitrary annotation dataframe
* remove chrY

# signals 0.11.4

* fix typo assigning A->B

# signals 0.11.3

* catch negative state_AS issue for singleton bins

# signals 0.11.2

* fix negative state_AS issue. This was caused during the fill missing step which
assigns states with NA BAF values based on neighbouring bins. When there was a 
state transition between bins sometimes A+B>state.

# signals 0.11.1

* for umap_clustering, default is now to not use PCA

# signals 0.11.0

* updates to plotting
* male/female option
* filter some cell for second pass phasing but do not remove
* do not remove cells by default if they have low coverage
* keep cells with large hom-dels

# signals 0.10.0

* add chr string check
* remove acrocentric chromosomes in arm consensus
* minor changes to heatmap and plotting

# signals 0.9.1

* Update docker to install suggests packages

# signals 0.9.0

* Fix phasing issue that happens when the cluster identified to phase relative to has a diploid region

# signals 0.8.0

* Add option to mask bins during inference

# signals 0.7.6

* Fix bug in plot_clusters_used_for_phasing

# signals 0.7.5

* Fixed r cmd check

# signals 0.7.4

# signals 0.7.3

# signals 0.7.2

* release for zenodo

# signals 0.7.1

# signals 0.7.0

* change name to signals
* multiple updates to plotting
* rewrite of scRNAseq

# signals 0.6.2

* fix missing argument

# signals 0.6.1

* changed colours
* fill in missing bins
* improved documentation about inputs

# signals 0.6.0

* add option to remove noisy cells from phasing

# signals 0.5.5

* update plotting and clustering

# signals 0.5.4

* update docker

# signals 0.5.3

* some small updates to plotting and default params

# signals 0.5.2

* update to vignette and plotting

# signals 0.5.1

* updated default parameters

# signals 0.5.0

* update ascn inference

# signals 0.4.3

* filtering utility function

# signals 0.4.2

* fix bug with filtering

# signals 0.4.1

* Fix bug in HMM when total copy number = 0
* Add function to filter hscn object

# signals 0.4.0

* Add option to filter haplotypes
* Added Dockerfile and github actions to push to Dockerhub
* Added QC metadata table to output

# signals 0.3.0

* Version associated with preprint
* some fixes to plotting

# signals 0.2.1

* updates to heatmap plotting

# signals 0.2.0

# signals 0.1.0

* Added a `NEWS.md` file to track changes to the package.
