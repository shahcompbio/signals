#' @export
assign_haplotype_label <- function(haplotypes, hapbinsize = 50e3){
  colnames(haplotypes) <- c("chr", "pos", "cell_id", "allele0", "allele1")
  haplotypes$chr <- as.character(haplotypes$chr)
  haplotypes$pos2 <- haplotypes$pos
  
  bins <- getBins(binsize = hapbinsize) %>%
    dplyr::rename(start_bins = start, end_bins = end, chr_bins = chr) %>%
    dplyr::select(-width) %>%
    as.data.table() %>%
    .[, hap_label := 1:.N]
  
  haplotypes <- haplotypes[bins, on = .(chr == chr_bins, pos > start_bins, pos < end_bins)]
  haplotypes <- na.omit(haplotypes)
  
  hap_labels <- dplyr::distinct(haplotypes, chr, pos2, hap_label) %>% dplyr::rename(position = pos2)
  return(hap_label)
}

#' @export
assign_label_persnp <- function(haplotypes, hapbinsize = 50e3){
  snpdf <- as.data.table(haplotypes)[, hap_label := .GRP, by = list(chr, position)]
  return(snpdf)
}

#' @export
assign_bins_haplotypes <- function(haplotypes, binsize = 0.5e6){
  binshaps <- haplotypes %>%
    .[, start := floor(position / binsize) * binsize + 1] %>%
    .[, end := start + binsize - 1] %>%
    .[, position := NULL]
  
  return(binshaps)
}
