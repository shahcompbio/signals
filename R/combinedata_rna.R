#' @export
phase_haplotypes_rna <- function(haplotypes) {
  phased_haplotypes <- data.table::as.data.table(haplotypes) %>%
    .[, lapply(.SD, sum), by = .(chr, hap_label), .SDcols = c("allele1", "allele0")] %>%
    .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
    .[, c("allele1", "allele0") := NULL]

  return(phased_haplotypes)
}

#' @export
format_haplotypes_rna <- function(haplotypes,
                                  filtern = 0,
                                  phased_haplotypes = NULL,
                                  phasing_method = "distribution",
                                  create_cell_id = FALSE,
                                  ...) {
  message("Phase haplotypes...")
  if (is.null(phased_haplotypes)) {
    if (phasing_method == "distribution") {
      message("Phasing based on distribution across all cells")
      phased_haplotypes <- phase_haplotypes_rna(haplotypes)
    }
    message("Join phased haplotypes...")
    haplotypes <- as.data.table(haplotypes)[phased_haplotypes,
      on = .(chr, hap_label)
    ] %>%
      .[!is.na(cell_id)]
  } else {
    message("Join phased haplotypes...")
    phased_haplotypes <- as.data.table(dplyr::distinct(phased_haplotypes, chr, hap_label, phase))
    haplotypes <- as.data.table(haplotypes)[phased_haplotypes,
      on = .(chr, hap_label), allow.cartesian = TRUE
    ] %>%
      .[!is.na(cell_id)]
  }

  message("Reorder haplotypes based on phase...")
  haplotypes <- haplotypes %>%
    .[, alleleA := ifelse(phase == "allele0", allele1, allele0)] %>%
    .[, alleleB := ifelse(phase == "allele0", allele0, allele1)] %>%
    .[, totalcounts := alleleA + alleleB] %>%
    .[totalcounts > filtern, BAF := alleleB / totalcounts] %>%
    .[, c("allele1", "allele0", "phase") := NULL]
  
  haplotype_snp_counts <-  haplotypes %>%
    .[, list(
      alleleA = sum(alleleA),
      alleleB = sum(alleleB)
    ), by = c("cell_id", "chr", "position", "hap_label")] %>%
    .[, totalcounts := alleleA + alleleB] %>%
    .[totalcounts > filtern, BAF := alleleB / totalcounts] %>%
    dplyr::rename(start = position) %>% 
    dplyr::select(cell_id, chr, start, hap_label, alleleA, alleleB, totalcounts, BAF)

  haplotypes <- haplotypes %>%
    .[, list(
      alleleA = sum(alleleA),
      alleleB = sum(alleleB),
      start = min(position), end = max(position)
    ), by = c("cell_id", "chr", "hap_label")] %>%
    .[, totalcounts := alleleA + alleleB] %>%
    .[totalcounts > filtern, BAF := alleleB / totalcounts] %>%
    dplyr::select(cell_id, chr, start, end, hap_label, alleleA, alleleB, totalcounts, BAF)

  if (create_cell_id) {
    haplotypes <- haplotypes %>% 
      dplyr::mutate(cell_id = paste(patient, sample, cell_id, sep = "-")) %>%
      as.data.frame() %>%
      orderdf()
    
    haplotype_snp_counts <- haplotype_snp_counts %>% 
      dplyr::mutate(cell_id = paste(patient, sample, cell_id, sep = "-")) %>%
      as.data.frame() %>%
      orderdf()
  }

  return(list(block_counts = haplotypes, snp_counts = haplotype_snp_counts))
}
