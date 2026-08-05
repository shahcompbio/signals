## Documentation for the datasets shipped in data/.
##
## Provenance is recorded where it could be established from how the data is used
## in the package. Several of the reference tables predate this documentation and
## their original source is not recorded anywhere in the repository; those have an
## empty @source rather than a guess.

# ---- example single cell data ------------------------------------------------

#' Copy number bins for DLP+ cells
#'
#' Total copy number per 500kb bin per cell, as produced by the HMMcopy step of
#' the DLP+ pipeline. Covers three OV2295-derived cell lines: SA921 (85 cells),
#' SA1090 (110) and SA922 (55).
#'
#' This is the format the copy number functions expect: `chr`, `start`, `end`,
#' `cell_id`, `state` (total copy number) and `copy` (normalised read depth).
#'
#' @format A data frame with 1,093,750 rows and 7 columns:
#' \describe{
#'   \item{chr, start, end}{genomic bin}
#'   \item{reads}{reads in the bin}
#'   \item{copy}{normalised read depth}
#'   \item{state}{total copy number}
#'   \item{cell_id}{cell identifier, `{sample}-{library}-R{row}-C{column}`}
#' }
#'
#' @seealso [callHaplotypeSpecificCN()], [plotCNprofile()], [SVs]
#' @md
"CNbins"

#' Haplotype block allele counts for DLP+ cells
#'
#' Raw (unphased) allele counts per haplotype block per cell for the same cells
#' as [CNbins]. Input to [format_haplotypes_dlp()] and then
#' [callHaplotypeSpecificCN()].
#'
#' @format A data frame with 13,268,159 rows and 7 columns:
#' \describe{
#'   \item{cell_id}{cell identifier}
#'   \item{chr, start, end}{haplotype block position}
#'   \item{allele_id}{0 or 1, distinguishing the two alleles}
#'   \item{hap_label}{haplotype block identifier, consistent across cells}
#'   \item{readcount}{reads supporting this allele in this cell}
#' }
#'
#' @seealso [format_haplotypes_dlp()], [callHaplotypeSpecificCN()]
#' @md
"haplotypes"

#' Per-cell quality control metrics
#'
#' HMMcopy and alignment QC metrics from the DLP+ pipeline, one row per cell.
#' `quality` is the classifier score used to filter cells, and `is_s_phase`
#' flags replicating cells.
#'
#' @format A data frame with 237 rows and 66 columns. Columns include
#'   `cell_id`, `quality`, `is_s_phase`, `total_mapped_reads`, `coverage_depth`,
#'   `mad_hmmcopy`, `breakpoints` and `state_mode`, plus the chip position
#'   (`row`, `column`), library and index metadata.
#' @md
"CNmetrics"

#' SNV counts for DLP+ cells
#'
#' Per-cell reference and alternate read counts at SNV positions, for the same
#' cells as [CNbins].
#'
#' @format A data frame with 585,531 rows and 9 columns:
#' \describe{
#'   \item{chr, start}{SNV position}
#'   \item{ref, alt}{reference and alternate base}
#'   \item{ref_counts, alt_counts, total_counts}{read counts in this cell}
#'   \item{cell_id, sample_id}{cell and sample identifiers}
#' }
#' @md
"SNV"

# ---- example scRNA data ------------------------------------------------------

#' Gene expression counts matrix
#'
#' Sparse gene-by-cell counts matrix for 10x scRNA-seq of the OV2295 lines, used
#' in the allele specific copy number from RNA vignette.
#'
#' @format A `dgCMatrix` with 21,175 genes (rows) and 814 cells (columns).
#'   Row names are gene symbols, column names are `{patient}-{sample}_{lane}_{barcode}`.
#' @md
"countsmatrix"

#' Haplotype allele counts for scRNA cells
#'
#' Per-cell allele counts at phased SNP positions for the scRNA-seq data,
#' matching [countsmatrix]. Input to the RNA copy number pipeline.
#'
#' @format A data frame with 1,049,775 rows and 10 columns:
#' \describe{
#'   \item{cell_id}{cell identifier}
#'   \item{chr, position}{SNP position}
#'   \item{ref, alt}{reference and alternate base}
#'   \item{hap_label}{haplotype block identifier}
#'   \item{allele0, allele1}{counts for each phased allele}
#'   \item{patient, sample}{sample metadata}
#' }
#' @md
"haplotypes_rna"

# ---- reference tables --------------------------------------------------------

#' Cytoband coordinates by genome build
#'
#' Cytoband (ideogram) coordinates, used to draw the ideogram track in
#' [plotCNprofile()] and [plotHeatmap()]. Indexed by genome build.
#'
#' @format A named list of 5 data tables (`hg38`, `hg19`, `hg18`, `hg17`,
#'   `hg16`), each with chromosome, start, end, band name and Giemsa stain
#'   level. The hg19 table has 1,293 rows.
#'
#' @source
#' @keywords internal
#' @md
"cytoband_map"

#' Gene coordinates
#'
#' Gene positions used by the `genes` argument of [plotCNprofile()] and
#' [plotHeatmap()] to annotate named genes.
#'
#' @format A named list with one element, `hg19`: a data frame of 23,657 rows
#'   with `chr`, `start`, `end` and `ensembl_gene_symbol`.
#'
#' @source
#' @keywords internal
#' @md
"gene_locations"

#' DLP+ 500kb bin coordinates
#'
#' The standard 500kb bin set for DLP+ data, used by the simulation functions to
#' generate copy number profiles on the same grid as real data.
#'
#' @format A data frame with 4,375 rows and 3 columns: `chr`, `start`, `end`.
#'
#' @keywords internal
#' @md
"dlpbins"

#' hg19 chromosome lengths
#'
#' Chromosome lengths in base pairs, the default used by [getBins()] to build
#' fixed-width bins.
#'
#' @format A named integer vector of length 25, names being chromosome names.
#'
#' @keywords internal
#' @md
"hg19_chrlength"

#' hg19 chromosome and arm coordinates
#'
#' Chromosome and chromosome-arm coordinates, used when mapping RNA data onto
#' arm-level segments.
#'
#' @format A data frame with 72 rows and 5 columns: `chr`, `start`, `end`, `arm`
#'   and `chrarm`.
#'
#' @source
#' @keywords internal
#' @md
"hg19chrom_coordinates"

#' Common and rare fragile sites, binned
#'
#' Fragile site annotation mapped onto 500kb bins, for annotating copy number
#' breakpoints against known fragile sites.
#'
#' @format A data frame with 6,260 rows and 6 columns:
#' \describe{
#'   \item{chr, start, end}{genomic bin}
#'   \item{frequency}{`common` (492 bins), `rare` (178) or `NA` (5,590)}
#'   \item{FSname}{fragile site name, where the bin overlaps one}
#'   \item{nFS}{number of fragile sites overlapping the bin}
#' }
#'
#' @source
#' @keywords internal
#' @md
"fragile_sites"

#' Fragile site locations
#'
#' The underlying fragile site intervals that [fragile_sites] bins.
#'
#' @format A data frame with 76 rows and 6 columns:
#' \describe{
#'   \item{seqnames, start, end}{fragile site interval, chromosomes `chr`-prefixed}
#'   \item{FSname}{fragile site name}
#'   \item{metadata}{description of the site}
#'   \item{frequency}{`common` or `rare`}
#' }
#'
#' @source
#' @keywords internal
#' @md
"fragile_sites_location"

# ---- structural variants -----------------------------------------------------

#' Structural variants for cell line SA921
#'
#' Structural variant breakpoints called by destruct on the OV2295-derived cell
#' line SA921, one of the samples in [CNbins], so the two can be plotted
#' together.
#'
#' Filtered to calls supported by at least 5 reads. destruct's own filters
#' (`is_filtered`, `is_germline`, `is_dgv`) were already applied upstream; the
#' read-support floor removes the long tail of low-confidence calls, most of
#' which are spurious foldbacks.
#'
#' This is the column format [plotCNprofile()] expects for its `SV` argument.
#' `strand_1` and `strand_2` must be `"+"`/`"-"`, and `read_count` is required
#' by the `"lines_and_arcs"` style.
#'
#' @format A data frame with 447 rows and 9 columns:
#' \describe{
#'   \item{chromosome_1, position_1, strand_1}{first breakend}
#'   \item{chromosome_2, position_2, strand_2}{second breakend}
#'   \item{type}{deletion, duplication, inversion or translocation}
#'   \item{rearrangement_type}{destruct's finer label, including foldback and balanced}
#'   \item{read_count}{number of supporting reads (destruct `num_reads`)}
#' }
#'
#' @source destruct breakpoint calls for SA921, from the OV2295 single cell
#'   whole genome sequencing dataset.
#'
#' @seealso [plotCNprofile()] for plotting these alongside copy number, and the
#'   "Structural variant visualization" vignette.
#' @md
"SVs"
