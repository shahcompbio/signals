#' Format haplotypes from DLP+ data
#'
#' Converts raw haplotype data from DLP+ sequencing to the format required by
#' signals. This includes binning haplotypes to match copy number bin coordinates
#' and converting from long to wide format.
#'
#' @param haplotypes A data.frame with raw haplotype allele counts. Required columns:
#'   `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele_id`, `readcount`.
#' @param CNbins A data.frame with copy number bin coordinates. Used to align haplotype
#'   bins. Required columns: `chr`, `start`, `end`.
#' @param hmmcopybinsize Bin size used by HMMcopy for copy number calling. Default 0.5e6 (500kb).
#'
#' @return A data.frame with formatted haplotypes containing columns:
#'   * `cell_id`: Cell identifier
#'   * `chr`: Chromosome
#'   * `start`, `end`: Bin coordinates (aligned to CNbins)
#'   * `hap_label`: Haplotype block identifier
#'   * `allele1`, `allele0`: Read counts for each allele
#'   * `totalcounts`: Total read counts (allele1 + allele0)
#'
#' @examples
#' data(CNbins)
#' data(haplotypes)
#' haps_formatted <- format_haplotypes_dlp(haplotypes, CNbins)
#'
#' @seealso [format_haplotypes()] for adding phasing information
#' @export
format_haplotypes_dlp <- function(haplotypes, CNbins, hmmcopybinsize = 0.5e6) {
  # Validate inputs
  raw_hap_cols <- c("cell_id", "chr", "start", "end", "hap_label", "allele_id", "readcount")
  check_required_columns(haplotypes, raw_hap_cols, "haplotypes")
  check_required_columns(CNbins, c("chr", "start", "end"), "CNbins")
  check_positive_numeric(hmmcopybinsize, "hmmcopybinsize")

  options("scipen" = 20)

  haplotypes <- haplotypes %>%
    data.table::as.data.table()

  CNbins <- CNbins %>%
    data.table::as.data.table()

  bins <- CNbins[, c("chr", "start", "end")] %>%
    unique(., by = c("chr", "start", "end")) %>%
    .[, binid := paste(chr, start, end, sep = "_")] %>%
    .$binid
  message(paste0("Number of distinct bins in copy number data: ", length(bins)))

  binshaps <- haplotypes[, c("chr", "start", "end")] %>%
    unique(., by = c("chr", "start", "end")) %>%
    .[, start := floor(start / hmmcopybinsize) * hmmcopybinsize + 1] %>%
    .[, end := start + hmmcopybinsize - 1] %>%
    .[, binid := paste(chr, start, end, sep = "_")] %>%
    .$binid

  message(paste0("Number of distinct bins in haplotype data: ", length(unique(binshaps))))

  formatted_haplotypes <- haplotypes %>%
    .[, allele_id := paste0("allele", allele_id)] %>%
    data.table::dcast(., ... ~ allele_id, value.var = "readcount", fill = 0L) %>%
    .[, start := floor(start / hmmcopybinsize) * hmmcopybinsize + 1] %>%
    .[, end := start + hmmcopybinsize - 1] %>%
    .[, lapply(.SD, sum), by = .(cell_id, chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
    .[, totalcounts := allele1 + allele0] %>%
    .[, hbinid := paste(chr, start, end, sep = "_")] %>%
    .[hbinid %in% bins] %>%
    .[, hbinid := NULL]

  binshaps2 <- formatted_haplotypes[, c("chr", "start", "end")] %>%
    unique(., by = c("chr", "start", "end")) %>%
    .[, start := floor(start / hmmcopybinsize) * hmmcopybinsize + 1] %>%
    .[, end := start + hmmcopybinsize - 1] %>%
    .[, binid := paste(chr, start, end, sep = "_")] %>%
    .$binid

  message(paste0("Number of distinct bins in formatted haplotype data: ", length(unique(binshaps2))))

  return(as.data.frame(formatted_haplotypes))
}

#' Format and phase haplotypes
#'
#' Adds phasing information to formatted haplotypes, determining which allele
#' is the "A" (dominant) allele and which is the "B" (minor) allele for each
#' haplotype block.
#'
#' @param haplotypes A data.frame with haplotype allele counts. Required columns:
#'   `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`, `totalcounts`.
#' @param filtern Minimum total read count per haplotype to include. Default 0.
#' @param hmmcopybinsize Bin size used by HMMcopy. Default 0.5e6 (500kb).
#' @param phased_haplotypes Optional pre-computed phasing from `computehaplotypecounts()`.
#'   If NULL, phasing is computed automatically using `phasing_method`.
#' @param phasing_method Method for phasing. "distribution" (default) uses the
#'   distribution across all cells. Alternative is to use top N imbalanced cells
#'   via `computehaplotypecounts()`.
#' @param ... Additional arguments passed to `computehaplotypecounts()` if
#'   `phasing_method` is not "distribution".
#'
#' @return A data.frame with phased haplotypes containing:
#'   * `alleleA`: Read counts for the dominant allele
#'   * `alleleB`: Read counts for the minor allele
#'   * `BAF`: B-allele frequency (alleleB / totalcounts)
#'
#' @seealso [format_haplotypes_dlp()] for initial formatting,
#'   [phase_haplotypes()] and [computehaplotypecounts()] for phasing methods
#' @export
format_haplotypes <- function(haplotypes,
                              filtern = 0,
                              hmmcopybinsize = 0.5e6,
                              phased_haplotypes = NULL,
                              phasing_method = "distribution", ...) {
  # Validate inputs
  validate_haplotypes(haplotypes, formatted = TRUE)
  check_positive_numeric(filtern, "filtern", allow_zero = TRUE)

  message("Phase haplotypes...")
  if (is.null(phased_haplotypes)) {
    if (phasing_method == "distribution") {
      message("Phasing based on distribution across all cells")
      phased_haplotypes <- phase_haplotypes(haplotypes)
    } else {
      phased_haplotypes <- computehaplotypecounts(haplotypes, ...)
    }
  }

  message("Join phased haplotypes...")
  haplotypes <- haplotypes[phased_haplotypes, on = .(chr, start, end, hap_label)]

  message("Reorder haplotypes based on phase...")
  haplotypes <- haplotypes %>%
    .[, alleleA := ifelse(phase == "allele0", allele1, allele0)] %>%
    .[, alleleB := ifelse(phase == "allele0", allele0, allele1)] %>%
    .[totalcounts > filtern, BAF := alleleB / totalcounts] %>%
    .[, c("allele1", "allele0", "phase") := NULL]

  return(haplotypes)
}

#' Phase haplotypes using distribution method
#'
#' Determines the phase (which allele is A vs B) for each haplotype block based on
#' the distribution of allele counts across all cells. The allele with fewer total
#' counts is assigned as the "A" allele.
#'
#' @param haplotypes A data.frame with formatted haplotype data. Required columns:
#'   `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`.
#'
#' @return A data.frame with phasing information containing:
#'   * `chr`, `start`, `end`, `hap_label`: Haplotype block identifiers
#'   * `phase`: Which allele ("allele0" or "allele1") should be treated as alleleA
#'
#' @seealso [computehaplotypecounts()] for alternative phasing using top imbalanced cells
#' @export
phase_haplotypes <- function(haplotypes) {
  phased_haplotypes <- data.table::as.data.table(haplotypes) %>%
    .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
    .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
    .[, c("allele1", "allele0") := NULL]

  return(phased_haplotypes)
}

#' Compute phasing using top imbalanced cells
#'
#' Determines the phase for each haplotype block using the top N cells with the
#' highest allelic imbalance. This method is more robust when there are many
#' balanced cells that could add noise to phasing.
#'
#' @param haplotypes A data.frame with formatted haplotype data. Required columns:
#'   `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`, `cell_id`.
#' @param ncells Number of top imbalanced cells to use for phasing. Default 10.
#' @param arm If TRUE, perform phasing per chromosome arm rather than whole chromosome.
#'   Default FALSE.
#'
#' @return A data.frame with phasing information containing:
#'   * `chr`, `start`, `end`, `hap_label`: Haplotype block identifiers
#'   * `phase`: Which allele should be treated as alleleA
#'
#' @seealso [phase_haplotypes()] for phasing using all cells
#' @export
computehaplotypecounts <- function(haplotypes, ncells = 10, arm = FALSE) {
  # Validate inputs
  check_positive_integer(ncells, "ncells")

  message(paste0("Phasing based on distribution across top ", ncells, " cells with highest imbalance"))
  formatted_haplotypes <- haplotypes %>%
    .[, R := fifelse(allele0 == 0 | allele1 == 0, 0, 1)] %>%
    .[, R0 := fifelse(allele0 == 0, "allele0", "allele1")] %>%
    .[, R1 := fifelse(allele1 == 0, "allele1", "allele0")]

  if (arm == FALSE) {
    perchr <- formatted_haplotypes %>%
      .[, list(meanR = mean(R), meanR0 = Mode(R0), meanR1 = Mode(R1)),
        by = c("cell_id", "chr")
      ] %>%
      # .[, dominant := fifelse(meanR0 < meanR1, "R1", "R0")] %>%
      setkey("meanR") %>%
      .[, head(.SD, ncells), by = c("chr")]
  } else {
    formatted_haplotypes$arm <- coord_to_arm(formatted_haplotypes$chr, formatted_haplotypes$start)
    perchr <- formatted_haplotypes %>%
      .[, list(meanR = mean(R), meanR0 = Mode(R0), meanR1 = Mode(R1)),
        by = c("cell_id", "chr", "arm")
      ] %>%
      # .[, dominant := fifelse(meanR0 < meanR1, "R1", "R0")] %>%
      setkey("meanR") %>%
      .[, head(.SD, ncells), by = c("chr", "arm")]
  }


  if (arm == FALSE) {
    limitedhaps <- perchr[formatted_haplotypes, on = c("cell_id", "chr"), nomatch = 0] %>%
      .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")]
    phased_haplotypes <- limitedhaps %>%
      dplyr::group_by(chr, start, end, hap_label) %>%
      dplyr::summarise(phase = Mode(phase), meanR = mean(meanR)) %>%
      dplyr::ungroup() %>%
      data.table::as.data.table()
  } else {
    limitedhaps <- perchr[formatted_haplotypes, on = c("cell_id", "chr", "arm"), nomatch = 0] %>%
      .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")]
    phased_haplotypes <- limitedhaps %>%
      dplyr::group_by(chr, start, arm, end, hap_label) %>%
      dplyr::summarise(phase = Mode(phase), meanR = mean(meanR)) %>%
      dplyr::ungroup() %>%
      data.table::as.data.table()
  }

  return(phased_haplotypes)
}

#' Combine copy number bins with haplotype data
#'
#' Joins copy number bin data with phased haplotype counts to produce a combined
#' data.frame with copy number states and B-allele frequency (BAF) for each bin.
#'
#' @param haplotypes A data.frame with haplotype allele counts. Required columns:
#'   `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele_id`, `readcount` (raw) or
#'   `allele1`, `allele0`, `totalcounts` (formatted).
#' @param CNbins A data.frame with copy number bin data. Required columns:
#'   `cell_id`, `chr`, `start`, `end`, `state`, `copy`.
#' @param filtern Minimum total read count per bin to include. Default 0.
#' @param phased_haplotypes Optional pre-computed phased haplotypes from
#'   `computehaplotypecounts()`. If NULL, phasing is performed automatically.
#' @param minbins Minimum number of bins per cell to include. Default 100.
#' @param minbinschr Minimum number of bins per chromosome per cell. Default 10.
#' @param phasing_method Method for phasing haplotypes. One of "distribution" (default)
#'   or use top N imbalanced cells.
#' @param ... Additional arguments passed to `format_haplotypes()`.
#'
#' @return A data.frame with columns from CNbins plus:
#'   * `alleleA`: Read counts for A allele
#'   * `alleleB`: Read counts for B allele
#'   * `totalcounts`: Total read counts
#'   * `BAF`: B-allele frequency (alleleB / totalcounts)
#'
#' @examples
#' data(CNbins)
#' data(haplotypes)
#' # Format haplotypes first
#' haps_formatted <- format_haplotypes_dlp(haplotypes, CNbins)
#' # Combine with CNbins
#' cnbaf <- combineBAFCN(haps_formatted, CNbins)
#'
#' @export
combineBAFCN <- function(haplotypes,
                         CNbins,
                         filtern = 0,
                         phased_haplotypes = NULL,
                         minbins = 100,
                         minbinschr = 10,
                         phasing_method = "distribution", ...) {
  # Validate inputs
  validate_cnbins(CNbins)
  check_positive_integer(minbins, "minbins", allow_zero = TRUE)
  check_positive_integer(minbinschr, "minbinschr", allow_zero = TRUE)
  check_positive_numeric(filtern, "filtern", allow_zero = TRUE)

  message("Finding overlapping cell IDs between CN data and haplotype data...")
  cellidoverlap <- intersect(CNbins$cell_id, haplotypes$cell_id)
  message(paste0("Total number of cells in both CN and haplotypes: ", length(cellidoverlap)))

  CNbins <- data.table::as.data.table(CNbins)
  haplotypes <- data.table::as.data.table(haplotypes)

  if (all(cellidoverlap %in% CNbins$cell_id)) {
    message(paste0("Number of cells in CN data: ", length(unique(CNbins$cell_id))))
    CNbins <- CNbins[cell_id %in% cellidoverlap]
  }

  if (all(cellidoverlap %in% haplotypes$cell_id)) {
    message(paste0("Number of cells in haplotype data: ", length(unique(CNbins$cell_id))))
    haplotypes <- haplotypes[cell_id %in% cellidoverlap]
  }

  message("Joining bins and haplotypes...")
  haplotypes <- format_haplotypes(haplotypes,
    phased_haplotypes = phased_haplotypes,
    phasing_method = phasing_method, ...
  )
  haplotypes <- data.table::as.data.table(haplotypes)

  CNbins <- CNbins[haplotypes, on = c("chr", "start", "end", "cell_id"), nomatch = 0]

  if (nrow(CNbins) == 0) {
    warning("Join between CNbins and haplotypes resulted in 0 rows. ",
            "Check that chr, start, end coordinates match between datasets.")
  }

  CNBAF <- data.table::as.data.table(CNbins) %>%
    .[totalcounts > filtern] %>%
    .[, lapply(.SD, sum), by = .(chr, start, end, cell_id, state, copy), .SDcols = c("alleleA", "alleleB", "totalcounts")] %>%
    .[, BAF := alleleB / totalcounts]

  CNBAF <- CNBAF %>%
    .[, n := .N, by = "cell_id"] %>%
    .[n > minbins] %>%
    .[, n := NULL] %>%
    .[, n := .N, by = c("chr", "cell_id")] %>%
    .[, n := ifelse(chr == "X", 100, n)] %>% 
    .[, ncell := min(n), by = c("cell_id")] %>%
    #do not filter out if chr X below cutoff
    .[ncell > minbinschr] %>%
    .[, n := NULL] %>%
    .[, ncell := NULL]

  message(paste0("Total number of cells after removing cells with < ", minbins, " bins: ", length(unique(CNBAF$cell_id))))

  return(as.data.frame(CNBAF))
}
