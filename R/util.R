Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' @export
createCNmatrix <- function(CNbins, field = "state", maxval = 11, na.rm = FALSE, fillna = FALSE){

  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))

  CNbins <- data.table::as.data.table(CNbins)

  cnmatrix <- CNbins %>%
    .[, segid := paste(chr, as.integer(start), as.integer(end), sep = "_")] %>%
    .[, state := data.table::fifelse(state > maxval, maxval, state)] %>%
    .[, width := end - start] %>%
    data.table::dcast(., chr + start + end + width ~ cell_id, value.var = field, fill = NA) %>%
    .[dfchr, on = "chr"] %>%
    .[order(idx, start)]

  if (fillna == TRUE){
  cnmatrix <- cnmatrix %>%
    dplyr::as_tibble() %>%
    tidyr::fill(., names(cnmatrix), .direction = "updown")
  }

  cnmatrix <- as.data.frame(cnmatrix)

  if (na.rm == TRUE){
    cnmatrix <- na.omit(cnmatrix)
  }

  rownames(cnmatrix) <- paste(cnmatrix$chr, as.integer(cnmatrix$start), as.integer(cnmatrix$end), sep = "_")
  cnmatrix = subset(cnmatrix, select = -c(idx))
  return(cnmatrix)
}

#' @export
createbreakpointmatrix <- function(segs, transpose = FALSE, internalonly = TRUE, use_state = FALSE){

  options("scipen"=20)

  if (use_state == FALSE){
    if (internalonly == TRUE){
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, row_num := .I] %>% # add row numbers
        .[, remove_row_num := .I[.N], by=.(cell_id, chr)]  %>% # find last row in each cell_id - chr group
        .[row_num != remove_row_num] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, tipInclusionProbabilities)
    }
  } else{
    if (internalonly == TRUE){
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, state, sep = "_")] %>%
        .[, row_num := .I] %>% # add row numbers
        .[, remove_row_num := .I[.N], by=.(cell_id, chr)]  %>% # find last row in each cell_id - chr group
        .[row_num != remove_row_num] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, state, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, tipInclusionProbabilities)
    }
  }

  segs_mat <- segs_bin %>%
    data.table::dcast(., loci ~ cell_id, value.var = "tipInclusionProbabilities", fill = 0)

  segs_mat <- as.data.frame(segs_mat)
  rownames(segs_mat) <- segs_mat$loci

  if (transpose == TRUE){
    segs_mat <- subset(segs_mat, select = -c(loci))
    segs_mat <- t(segs_mat)
  }

  return(segs_mat)
}


#' Make fixed-width bins
#'
#' Make fixed-width bins based on given bin size.
#'
#' @param chrom.lengths A named character vector with chromosome lengths. Names correspond to chromosomes.
#' @param binsize Size of bins in basepairs
#' @param chromosomes A subset of chromosomes for which the bins are generated.
#' @return A data frame with chr start and end positions
#' @export
getBins <- function(chrom.lengths=hg19_chrlength, binsize=1e6, chromosomes=NULL) {

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package \"GenomeInfoDb\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
  chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  ## Stop if none of the specified chromosomes exist
  if (length(chroms2use)==0) {
    chrstring <- paste0(chromosomes, collapse=', ')
    stop('Could not find length information for any of the specified chromosomes: ', chrstring)
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff)>0) {
    diffs <- paste0(diff, collapse=', ')
    warning('Could not find length information for the following chromosomes: ', diffs)
  }

  ### Making fixed-width bins ###

  message("Making fixed-width bins for bin size ", binsize, " ...")
  chrom.lengths.floor <- floor(chrom.lengths / binsize) * binsize
  clfloor2use <- chrom.lengths.floor[chroms2use]
  clfloor2use <- clfloor2use[clfloor2use >= binsize]
  if (length(clfloor2use) == 0) {
    stop("All selected chromosomes are smaller than binsize ", binsize)
  }
  bins <- unlist(GenomicRanges::tileGenome(clfloor2use, tilewidth=binsize), use.names=FALSE)

  GenomeInfoDb::seqlevels(bins) <- chroms2use
  GenomeInfoDb::seqlengths(bins) <- chrom.lengths[GenomeInfoDb::seqlevels(bins)]
  skipped.chroms <- setdiff(chromosomes, as.character(unique(GenomeInfoDb::seqnames(bins))))
  bins <- GenomeInfoDb::dropSeqlevels(bins, skipped.chroms, pruning.mode = 'coarse')

  if (length(skipped.chroms)>0) {
    warning("The following chromosomes are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse=', '))
  }


  bins <- as.data.frame(bins) %>%
    dplyr::select(-strand) %>% dplyr::rename(chr = seqnames)

  if (any(bins$width!=binsize)) {
    stop("tileGenome failed")
  }

  return(bins)
}

#' @export
widen_haplotypebins <- function(haplotypes, binsize = 5e6){
  message("Create GRanges objects...")
  bins <- getBins(binsize = binsize)
  bins_g <- GenomicRanges::makeGRangesFromDataFrame(bins, keep.extra.columns = TRUE)
  haplotypes_g <- GenomicRanges::makeGRangesFromDataFrame(haplotypes, keep.extra.columns = TRUE)

  message("Find overlaps...")
  overlaps <- GenomicRanges::findOverlaps(haplotypes_g, bins_g, ignore.strand = TRUE)

  message("Create new bins coordinates dataframe")
  bincoords <- data.table::as.data.table(bins_g[S4Vectors::subjectHits(overlaps)]) %>%
    .[, strand := NULL]

  message("Create haplotypes dataframe and merge with new bin coordinates...")
  newhaps <- data.table::as.data.table(haplotypes_g[S4Vectors::queryHits(overlaps)]) %>%
    .[, c("strand", "start", "seqnames", "end", "width") := NULL] %>%
    dplyr::bind_cols(., bincoords) %>%
    data.table::setnames(., "seqnames", "chr") %>%
    data.table::setcolorder(., c("chr", "start", "end", "cell_id"))

  return(newhaps)
}

#' @export
widen_bins <- function(CNbins,
                       binsize = 10e6,
                       roundstate = TRUE,
                       genome_coords = hg19_chrlength){
  message("Create GRanges objects...")
  bins <- getBins(genome_coords, binsize = binsize)
  bins_g <- GenomicRanges::makeGRangesFromDataFrame(bins, keep.extra.columns = TRUE)
  CNbinstemp_g <- GenomicRanges::makeGRangesFromDataFrame(as.data.frame(CNbins), keep.extra.columns = TRUE)

  message("Find overlaps...")
  overlaps <- GenomicRanges::findOverlaps(CNbinstemp_g, bins_g, ignore.strand = TRUE)

  message("Create new bins coordinates dataframe")
  bincoords <- data.table::as.data.table(bins_g[S4Vectors::subjectHits(overlaps)]) %>%
    .[, strand := NULL]

  message("Create CN dataframe and merge with new bin coordinates...")
  widerCNbins <- data.table::as.data.table(CNbinstemp_g[S4Vectors::queryHits(overlaps)]) %>%
    .[, c("strand", "start", "seqnames", "end", "width") := NULL] %>%
    dplyr::bind_cols(., bincoords) %>%
    data.table::setnames(., "seqnames", "chr") %>%
    data.table::setcolorder(., c("chr", "start", "end", "cell_id")) %>%
    .[, .(state = round(mean(state, na.rm = TRUE)),
          copy = mean(copy, na.rm = TRUE)), by = .(chr, start, end, cell_id)]

  return(as.data.frame(widerCNbins))
}

#' @export
snv_states <- function(SNV, CNbins){

  CN <- CNbins %>%
    dplyr::rename(chry = chr, starty = start, endy = end, cell_idy = cell_id) %>%
    as.data.table()

  SNV <- as.data.table(SNV)

  mappedSNVs <- CN[SNV,
                   on = .(chry == chr, cell_idy == cell_id, starty < start, endy > start)
                   ]
  mappedSNVs <- mappedSNVs %>%
    .[, end := NULL] %>%
    data.table::setnames(., "chry", "chr") %>%
    data.table::setnames(., "starty", "start") %>%
    data.table::setnames(., "cell_idy", "cell_id") %>%
    data.table::setcolorder(., c("chr", "start","ref", "alt", "cell_id")) %>%
    .[order(cell_id, chr, start)]

  return(as.data.frame(mappedSNVs))
}


#' Genomic coordinate to chromosome arm
#'
#' Returns chromosome arms for given chromosome and genomic position.
#'
#' @param chromosome Character or numeric vector, with chromosome of genomic coordinate
#' @param position Numeric vector, with genomic position within chromosome
#' @param assembly a string specifying which genome assembly version should be applied
#'   to determine chromosome arms. Allowed options are "hg38", hg19", "hg18", "hg17"
#'   and "hg16" (corresponding to the five latest human genome annotations in the
#'   UCSC genome browser).
#' @return Character vector, with choromosome arm of given genomic coordinates
coord_to_arm <- function(chromosome, position, assembly = "hg19", full = F){
  if(length(chromosome) !=  length(position)){
    stop("chromosome and position must have equal length")
  }
  if (!(assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16"))) {
    stop("Invalid assembly, allowed options are hg38, hg19, hg18, hg17 and hg16")
  }
  if(any(substr(chromosome, 1, 3) != "chr")){
    chromosome <- paste0("chr", chromosome)
  }
  if(any(!grepl("chr[X-Y]|[0-9]+", chromosome))){
    stop("Invalid chromosome, must be 1-22, X or Y (or chr1-chr22, chrX or chrY)")
  }
  data("cytoband_map", envir = environment())
  arms <- rep("     ", length(chromosome))
  for(i in unique(chromosome)){
    map <- cytoband_map[[assembly]][V1 == i]
    arm <- map[(findInterval(position[chromosome == i], map$V3)+1)]$V4
    if(!full){
      arm <- substr(arm, 1,1)
    }
    arms[chromosome == i] <- arm
  }
  return(arms)
}

#' @export
create_segments <- function(CNbins, field = "state"){
  newsegs <- CNbins %>%
    data.table::as.data.table() %>%
    .[order(cell_id, chr, start)] %>%
    .[, rlid := data.table::rleid(get(field)), by = cell_id] %>%
    .[, list(start = min(start),
             end = max(end)), by = .(cell_id, chr, get(field), rlid)] %>%
    .[order(cell_id, chr, start)] %>%
    dplyr::select(cell_id, chr, start, end, dplyr::everything(), -rlid) %>%
    as.data.frame()
  setnames(newsegs, "get", field)
  return(newsegs)
}

#' @export
orderdf <- function(CNbins){
  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))
  return(CNbins %>%
    as.data.table() %>%
    .[dfchr, on = "chr"] %>%
    .[order(cell_id, idx, start)] %>%
    .[, idx := NULL] %>%
      as.data.frame())
}

densmode <- function(x){
  dens <- density(x)
  dens$x[which.max(dens$y)]
}

#' @export
qc_summary <- function(cn){
  if (is.hscn(cn) | is.ascn(cn)){
    cn <- cn$data
  } else{
    cn <- cn
  }

  distance_df <- cn %>%
    dplyr::filter(state_AS_phased != "0|0") %>%
    dplyr::group_by(state_AS_phased, Min, Maj) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frac = n / dplyr::n()) %>%
    dplyr::filter(n > 10) %>%
    dplyr::group_by(state_AS_phased, Min, Maj, n, frac) %>%
    dplyr::summarise(medianBAF = median(BAF, na.rm = TRUE),
                     meanBAF = mean(BAF, na.rm = TRUE),
                     #modeBAF = densmode(BAF),
                     high95 = quantile(BAF, 0.975, na.rm = TRUE),
                     low95 = quantile(BAF, 0.025, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(expBAF = Min / (Min + Maj)) %>%
    dplyr::mutate(distance = sqrt((expBAF - medianBAF)^2)) %>%
    as.data.frame()

  summary_distance <- weighted.mean(distance_df$distance, distance_df$frac)

  message(paste0("Average distance from median to expected BAF = ", round(summary_distance, 4)))

  return(list(distance = distance_df, summary = summary_distance))
}
