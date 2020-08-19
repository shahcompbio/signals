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
fixjitter <- function(bps, nextend = 2){
  x <- as.data.frame(colSums(bps))
  names(x) <- "frequency"
  x$loci <- row.names(x)

  frequency_order <- order(x$frequency, decreasing = TRUE)
  nloci <- dim(bps)[2]

  while(length(frequency_order) > 0){
    idx <- frequency_order[1]
    minidx <- max(1, idx - nextend)
    maxidx <- min(nloci, idx + nextend)
    temp_mat <-  bps[, minidx:maxidx] #extract matrix nextend either side of locus of interest
    bps[, minidx:maxidx] <- 0 #set matrix to 0
    bps[, idx] <- as.numeric(rowSums(temp_mat) > 0)
    frequency_order <- setdiff(frequency_order, minidx:maxidx)
  }

  message(paste0("Original number of loci: ", nloci))
  x <- colSums(bps)
  x <- x[x>0]
  bps <- bps[,names(x)]
  message(paste0("New number of loci: ", dim(bps)[2]))

  return(as.data.frame(bps))
}

#' @export
createbreakpointmatrix <- function(segs,
                                   transpose = FALSE,
                                   internalonly = TRUE,
                                   use_state = FALSE,
                                   state_remove = 2,
                                   fixjitter = FALSE){

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
        dplyr::select(cell_id, loci, chr, start, end, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[state != state_remove] %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, chr, start, end, tipInclusionProbabilities)
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
        dplyr::select(cell_id, loci, chr, start, end, tipInclusionProbabilities)
    } else {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[state != state_remove] %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, state, sep = "_")] %>%
        .[, tipInclusionProbabilities := 1] %>%
        dplyr::select(cell_id, loci, chr, start, end, tipInclusionProbabilities)
    }
  }

  mapping <- dplyr::select(segs_bin, cell_id, chr, start, end, loci)

  segs_mat <- segs_bin %>%
    dplyr::select(cell_id, loci, tipInclusionProbabilities) %>%
    data.table::dcast(., loci ~ cell_id, value.var = "tipInclusionProbabilities", fill = 0) %>%
    .[gtools::mixedorder(loci)]

  segs_mat <- as.data.frame(segs_mat)
  rownames(segs_mat) <- segs_mat$loci

  if (transpose == TRUE){
    segs_mat <- subset(segs_mat, select = -c(loci))
    segs_mat <- t(segs_mat)
    if (fixjitter == TRUE){
      segs_mat <- fixjitter(segs_mat, nextend = 2)
    }
  }

  return(list(bps = segs_mat, mapping = mapping))
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
coord_to_arm <- function(chromosome, position, assembly = "hg19", full = FALSE, mergesmallarms = FALSE){
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

  small <- cytoband_map[[assembly]] %>%
    dplyr::mutate(arm = substr(V4, 1,1)) %>%
    dplyr::group_by(V1, arm) %>%
    dplyr::summarise(start = min(V2), end = max(V3)) %>% dplyr::mutate(width = (end- start) /1e6) %>%
    dplyr::mutate(small = ifelse(width < 40, TRUE, FALSE)) %>%
    dplyr::group_by(V1) %>%
    dplyr::summarise(small = any(small))

  for(i in unique(chromosome)){
    map <- cytoband_map[[assembly]][V1 == i]
    arm <- map[(findInterval(position[chromosome == i], map$V3)+1)]$V4
    if(!full){
      if (mergesmallarms & small[small$V1 == i,]$small){
        arm <- ""
      } else{
        arm <- substr(arm, 1,1)
      }
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

#' @export
per_arm_baf_mat <- function(haps){
  baf <- haps %>%
    dplyr::filter(chr != "Y") %>%
    dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = TRUE)) %>%
    dplyr::mutate(chrarm = paste0(chr, arm)) %>%
    dplyr::group_by(chr, arm, chrarm, cell_id) %>%
    dplyr::summarise(alleleA = sum(alleleA, na.rm = TRUE),
                     alleleB = sum(alleleB, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BAF = alleleB / (alleleA + alleleB),
                  total = alleleB + alleleA) %>%
    dplyr::mutate(idx = ifelse(chr == "X", 2 * 23, 2 * as.numeric(chr))) %>%
    dplyr::mutate(idx = ifelse(arm == "p", idx - 1, idx))

  idx <- dplyr::distinct(baf, chr, arm, chrarm, idx) %>% dplyr::arrange(idx)

  baf_mat <- baf %>%
    dplyr::select(cell_id, chrarm, BAF) %>%
    tidyr::pivot_wider(names_from = "chrarm", values_from = "BAF") %>%
    as.data.frame()

  row.names(baf_mat) <- baf_mat$cell_id
  baf_mat = subset(baf_mat, select = -c(cell_id))

  baf_mat <- baf_mat[, idx$chrarm]

  return(list(bafperchr = baf, bafperchrmat = baf_mat))
}


#' @export
per_arm_cn <- function(hscn){
  hscn_arm <- hscn %>%
    dplyr::filter(chr != "Y") %>%
    dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = TRUE)) %>%
    dplyr::mutate(chrarm = paste0(chr, arm)) %>%
    dplyr::group_by(chr, arm, chrarm, cell_id) %>%
    dplyr::summarise(alleleA = sum(alleleA, na.rm = TRUE),
                     alleleB = sum(alleleB, na.rm = TRUE),
                     state = Mode(state),
                     copy = median(copy),
                     state_AS_phased = Mode(state_AS_phased),
                     state_AS = Mode(state_AS),
                     LOH = Mode(LOH),
                     state_BAF = Mode(state_BAF),
                     phase = Mode(phase)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(BAF = alleleB / (alleleA + alleleB),
                  total = alleleB + alleleA) %>%
    dplyr::mutate(idx = ifelse(chr == "X", 2 * 23, 2 * as.numeric(chr))) %>%
    dplyr::mutate(idx = ifelse(arm == "p", idx - 1, idx)) %>%
    dplyr::mutate(start = data.table::fifelse(arm == "p", 1, 10)) %>%
    dplyr::mutate(end = data.table::fifelse(arm == "p", 2, 11))

  return(hscn_arm)
}
