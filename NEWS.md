# signals 0.18.0

* Add a `seed` argument to `callHaplotypeSpecificCN` and `callAlleleSpecificCN`. Phasing has three stochastic steps — the subsampling in `min_cells` that sets the cluster size, the UMAP embedding used to choose which cells phase each chromosome, and the subsampling in the beta-binomial fit. Left unseeded, repeated runs on identical input can select different cells to phase a chromosome with and so return different haplotype-specific states. `seed` is threaded through `proportion_imbalance`, `get_cells_per_chr_local`, `get_cells_per_chr_global`, `min_cells` and `fitBB`; the default remains `NULL` (unseeded), so existing behaviour is unchanged.
* Pass `n_sgd_threads = 0` to `uwot::umap` in `umap_clustering` and `umap_clustering_breakpoints`. uwot's SGD is only reproducible single-threaded, so a seed alone does not pin the embedding if that default ever changes.
* Fix two `@param` names that did not match their arguments (`viterbver` -> `viterbiver`, `global_phasing_for_diploid` -> `global_phasing_for_balanced`), which left both arguments undocumented.
* Fix `phasebyarm = TRUE` producing an empty phasing table, which made arm-level phasing unusable on either selection path. `phase_haplotypes_bychr()` filters haplotypes on `chrarm`, so the per-unit cell list must be keyed by chromosome arm, but `proportion_imbalance()` accepted `phasebyarm` without forwarding it, and `get_cells_per_chr_local()` — the default, since `cluster_per_chr = TRUE` — accepted it and ignored it, keying the list by chromosome. Since `"4p" == "4"` is never true, zero blocks were phased and `format_haplotypes()` then dropped every haplotype block. `get_cells_per_chr_local()` now clusters and keys by the phasing unit, `proportion_imbalance()` forwards the argument, and `plot_clusters_used_for_phasing()` strips the arm suffix before filtering on `chr`. This matters where imbalance is confined to one arm: on a 725-cell test sample chr4 has 17 cells with p-arm LOH that phase the p arm correctly but are balanced on the q arm, so q-arm blocks were phased off noise and the one cell with a whole-chromosome gain read `2|1` over 94 of 100 p-arm bins while alternating `1|2`/`2|1` across 12 runs on the q arm.
* Fix `get_cells_per_chr_local()` — the default path, since `cluster_per_chr = TRUE` — ranking hdbscan's noise cluster alongside real clusters when choosing which cells phase a unit. `get_cells_per_chr_global()` already excluded it. Since a grab-bag of outliers can easily look the most imbalanced, the default path could phase a chromosome off exactly the cells least likely to share a copy number state. Cluster `"0"` is now dropped before ranking, unless every cell landed in it, in which case they are all kept.
* Fix total copy number being a function of the phasing. Under `fillmissing = TRUE`, a bin left with `A + B > state` had its `state` set to `NA` and refilled from a neighbour, so `state` in the output was not always the `state` supplied in `CNbins`. The alleles are now clamped instead (`B` into `[0, state]`, then `A = state - B`), leaving `state` exactly as given. On a 369-cell sample this takes bins where `A + B != state` from 269 to zero.
* Add `min_propA` to `callHaplotypeSpecificCN`, default `0` (off, existing behaviour). `get_cells_per_chr_local()` selects the cluster with the highest proportion of imbalanced bins however low that proportion is, so a chromosome where nothing is measurably imbalanced still received a confident phasing driven by noise in one small cluster. Above zero, such a unit pools every cell instead. The best proportion per unit is now recorded either way, so unphaseable units are identifiable rather than indistinguishable from well-phased ones.
* Fix blocks with no reads in the selected phasing cells being dropped from every cell in the sample. `phase_haplotypes_bychr()` pooled counts over the selected subset only, so such a block produced no row at all, and the right join in `format_haplotypes()` then discarded those bins cohort-wide rather than just for the cells used to phase — the smaller the phasing subset, the more bins were lost. Those blocks now fall back to the all-cell majority phase via the new internal `phase_with_fallback()`, used by all three phasing branches.
* Fix one haplotype block being given two different phases. A block is identified by `chr` + `hap_label`, but `phase_haplotypes()` and `phase_haplotypes_bychr()` pooled counts per (block, bin), and the 500kb binning splits a block across two bins whenever it straddles a boundary, so the two halves were phased independently and could disagree — 268 of 6,477 blocks spanned two bins on one test chromosome and 84 of those received conflicting phases. Counts are now pooled per block, the phase decided once, and mapped back onto every bin the block spans. Where `global_phasing_for_balanced = TRUE` splits a chromosome into balanced and imbalanced bins, each block is assigned to one side by a majority of its bins rather than being split across the two cell sets.
* Add `allele_aware` to `callHaplotypeSpecificCN`, default `FALSE` (the existing flat transition matrix). The HMM state is the B allele copy number with `A = state - B`, so a flat matrix over B makes "B unchanged" the cheap self-transition and gives "A unchanged" no standing: `1|1 -> 2|1` was cheap while `1|1 -> 1|2` was penalised, for no reason beyond how the state space is parameterised. Setting this scores a transition by how many alleles change, `1[dA != 0] + 1[dB != 0]`, which is symmetric in A and B. The cost depends on the change in total copy number between the two bins, so the matrix is position-dependent, via the new `viterbi_pd()`; wherever total copy number is constant it is identical to the flat matrix, so only transitions across a copy number change are affected, where a one-allele move is preferred - `2|1 -> 0|1` over `2|1 -> 1|0`. On two DLP+ samples this left parent-of-origin correctness unchanged (1748/1748 and 1630/1631) and reduced allele-specific segments by 1.3-1.5%.

# signals 0.17.0

* **Breaking:** fix the Viterbi backtrace in `viterbi()` (C++) and `viterbiR()`. It seeded the final bin with the predecessor of the best final state and then backtracked from the per-column argmax instead of following the stored backpointers, so it did not return the most likely path and emitted spurious single-bin state changes. Single-bin sequences also decoded to state 0 (C++) or errored (R). This changes haplotype- and allele-specific copy number calls: on a 725-cell DLP+ sample (chr6, analysis from #79) 194 cells changed path, 0.30% of bins changed, and mean A/B segments per cell fell from 1.71 to 0.26. Use signals 0.16.0 to reproduce earlier results.
* Fix `alleleHMM` producing all `-Inf` emissions for homozygous-deletion bins (total CN 0), which forced the rest of the chromosome to minor CN 0.
* Fix `createCNmatrix(centromere = TRUE)`: masking is now gated on `fillnaplot`, so it no longer errors with `fillna = TRUE` and actually runs with `fillnaplot = TRUE`.

# signals 0.16.0

* Add sv_arcs_above to plotCNprofile: draw SV arcs in a band above the copy number panel rather than on top of it, split by copy number effect with apex scaled by genomic span
* Fix SVs with both breakends in the same bin being silently dropped, which lost every foldback at 10kb resolution
* Fix panels with no SVs not reserving the SV band, so stacked panels stay aligned
* Add show_chrbreaks and ybreaks options to plotCNprofile
* Fix multi-region tick marks so each region gets a tick at its start, and dividers bracket both region edges
* Add multi-region plotting support to plotCNprofile
* Add y-axis transform option for the mean/IQR heatmap track
* Fix ordered_cell_ids initialisation in plotHeatmap
* Add SVs example data (destruct breakpoints for SA921, the sample CNbins comes from) and a structural variant visualization vignette
* Document all bundled datasets, which had no help pages

# signals 0.15.0

* Add configurable annotation colour overrides for discrete and continuous annotation columns
* Add configurable continuous/discrete threshold handling and robust cell matching for annotation metadata
* Add a mean + IQR top summary track, with support for sourcing that track from a different column than the heatmap
* Refine SV orientation colours and document the new plotHeatmap options

# signals 0.14.1

* Fix rephasebins stability to ensure stable results

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
