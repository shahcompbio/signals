#' Structural variants for cell line SA921
#'
#' Structural variant breakpoints called by destruct on the OV2295-derived cell
#' line SA921, the same sample the bundled [CNbins] copy number data comes from,
#' so the two can be plotted together.
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
"SVs"
