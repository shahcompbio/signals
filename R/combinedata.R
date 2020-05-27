#' @export
format_haplotypes_dlp <- function(haplotypes, CNbins, hmmcopybinsize = 0.5e6){

  bins <- dplyr::distinct(CNbins, chr, start, end) %>%
    dplyr::mutate(binid = paste(chr, start, end, sep = "_")) %>%
    dplyr::pull(binid)

  formatted_haplotypes <- haplotypes %>%
    data.table::as.data.table() %>%
    .[, allele_id := paste0("allele", allele_id)] %>%
    data.table::dcast(., ... ~ allele_id, value.var = "readcount", fill = 0L) %>%
    .[, start := round(start / hmmcopybinsize) * hmmcopybinsize + 1] %>%
    .[, end := start + hmmcopybinsize - 1] %>%
    .[, lapply(.SD, sum), by = .(cell_id, chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
    .[, totalcounts := allele1 + allele0] %>%
    .[, hbinid := paste(chr, start, end, sep = "_")] %>%
    .[hbinid %in% bins] %>%
    .[, hbinid := NULL]

  return(as.data.frame(formatted_haplotypes))
}

#' @export
format_haplotypes <- function(haplotypes,
                              filtern = 0,
                              hmmcopybinsize = 0.5e6,
                              phased_haplotypes = NULL,
                              phasing_method = "distribution", ...){
  message("Pivot data frame...")

  message("Phase haplotypes...")
  if (is.null(phased_haplotypes)){
    if (phasing_method == "distribution"){
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

#' @export
phase_haplotypes <- function(haplotypes){
  phased_haplotypes <- data.table::as.data.table(haplotypes) %>%
    .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
    .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
    .[, c("allele1", "allele0") := NULL]

  return(phased_haplotypes)
}

#' @export
computehaplotypecounts <- function(haplotypes, ncells = 10, arm = FALSE){
  message(paste0("Phasing based on distribution across top ", ncells," cells with highest imbalance"))
  formatted_haplotypes <- haplotypes %>%
    .[, R := fifelse(allele0 == 0 | allele1 == 0, 0, 1)] %>%
    .[, R0 := fifelse(allele0 == 0, "allele0", "allele1")] %>%
    .[, R1 := fifelse(allele1 == 0, "allele1", "allele0")]

  if (arm == FALSE){
    perchr <- formatted_haplotypes %>%
      .[, list(meanR = mean(R), meanR0 = Mode(R0), meanR1 = Mode(R1)),
        by = c("cell_id", "chr")] %>%
      #.[, dominant := fifelse(meanR0 < meanR1, "R1", "R0")] %>%
      setkey("meanR") %>%
      .[, head(.SD, ncells), by = c("chr")]
  } else{
    formatted_haplotypes$arm <- coord_to_arm(formatted_haplotypes$chr, formatted_haplotypes$start)
    perchr <- formatted_haplotypes %>%
      .[, list(meanR = mean(R), meanR0 = Mode(R0), meanR1 = Mode(R1)),
        by = c("cell_id", "chr", "arm")] %>%
      #.[, dominant := fifelse(meanR0 < meanR1, "R1", "R0")] %>%
      setkey("meanR") %>%
      .[, head(.SD, ncells), by = c("chr", "arm")]
  }


  if (arm == FALSE){
    limitedhaps <- perchr[formatted_haplotypes, on = c("cell_id", "chr"), nomatch = 0] %>%
      .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")]
    phased_haplotypes <- limitedhaps %>%
      dplyr::group_by(chr, start, end, hap_label) %>%
      dplyr::summarise(phase = Mode(phase), meanR = mean(meanR)) %>%
      dplyr::ungroup() %>%
      data.table::as.data.table()
  } else{
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

#' @export
combineBAFCN <- function(haplotypes,
                         CNbins,
                         filtern = 0,
                         phased_haplotypes = NULL,
                         minbins = 100,
                         minbinschr = 10,
                         phasing_method = "distribution", ...){

  message("Finding overlapping cell IDs between CN data and haplotype data...")
  cellidoverlap <- intersect(CNbins$cell_id, haplotypes$cell_id)
  message(paste0("Total number of cells in both CN and haplotypes: " , length(cellidoverlap)))

  CNbins <- data.table::as.data.table(CNbins)
  haplotypes <- data.table::as.data.table(haplotypes)

  if (all(cellidoverlap %in% CNbins$cell_id)){
    message(paste0("Number of cells in CN data: ", length(unique(CNbins$cell_id))))
    message("Removing cells from CN data...")
    CNbins <- CNbins[cell_id %in% cellidoverlap]
  }

  if (all(cellidoverlap %in% haplotypes$cell_id)){
    message(paste0("Number of cells in haplotype data: ", length(unique(CNbins$cell_id))))
    message("Removing cells from haplotype data...")
    haplotypes <- haplotypes[cell_id %in% cellidoverlap]
  }

  message("Reformatting haplotypes")
  haplotypes <- format_haplotypes(haplotypes,
                                  phased_haplotypes = phased_haplotypes,
                                  phasing_method = phasing_method, ...)
  haplotypes <- data.table::as.data.table(haplotypes)

  message("Joining bins and haplotypes...")
  #CNbins <- data.table::merge.data.table(CNbins, haplotypes)
  CNbins <- CNbins[haplotypes, on = c("chr", "start", "end", "cell_id"), nomatch=0]

  message("Calculate BAF per bin...")
  CNBAF <- data.table::as.data.table(CNbins) %>%
    .[totalcounts > filtern] %>%
    .[, lapply(.SD, sum), by = .(chr, start, end, cell_id, state, copy), .SDcols = c("alleleA", "alleleB", "totalcounts")] %>%
    .[, BAF := alleleB / totalcounts]

  CNBAF <- CNBAF %>%
    .[, n := .N, by = "cell_id"] %>%
    .[n > minbins] %>%
    .[, n := NULL] %>%
    .[, n := .N, by = c("chr", "cell_id")] %>%
    .[, ncell := min(n), by = c("cell_id")] %>%
    .[ncell > minbinschr] %>%
    .[, n := NULL] %>%
    .[, ncell := NULL]

  message(paste0("Total number of cells after removing cells with < " , minbins, " bins: ", length(unique(CNBAF$cell_id))))

  return(as.data.frame(CNBAF))
}
