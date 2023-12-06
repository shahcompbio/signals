Mode <- function(x) {
  ux <- unique(na.exclude(x))
  ux[which.max(tabulate(match(x, ux)))]
}

#' @export
createCNmatrix <- function(CNbins,
                           field = "state",
                           maxval = 11,
                           na.rm = FALSE,
                           fillna = FALSE,
                           fillnaplot = FALSE,
                           wholegenome = FALSE,
                           genome = "hg19",
                           centromere = FALSE) {
  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))
  dfchr <- dplyr::filter(dfchr, chr %in% unique(CNbins$chr))

  CNbins <- data.table::as.data.table(CNbins)

  if (wholegenome == TRUE){
    genome <- as.data.table(getBins(binsize = CNbins$end[1] - CNbins$start[1]+1))
    CNbins <- lapply(unique(CNbins$cell_id),
                function(x) merge(genome,
                                  CNbins[cell_id == x],
                                  on = c("chr", "start", "end"),
                                  all = TRUE)[, cell_id := x]) %>%
              rbindlist(.)
  }

  if (field == "state") {
    cnmatrix <- CNbins %>%
      .[, segid := paste(chr, as.integer(start), as.integer(end), sep = "_")] %>%
      .[, state := data.table::fifelse(state > maxval, maxval, state)] %>%
      .[, width := end - start] %>%
      data.table::dcast(., chr + start + end + width ~ cell_id, value.var = field, fill = NA) %>%
      .[dfchr, on = "chr"] %>%
      .[order(idx, start)]
  } else {
    cnmatrix <- CNbins %>%
      .[, segid := paste(chr, as.integer(start), as.integer(end), sep = "_")] %>%
      .[, width := end - start] %>%
      data.table::dcast(., chr + start + end + width ~ cell_id, value.var = field, fill = NA) %>%
      .[dfchr, on = "chr"] %>%
      .[order(idx, start)]
  }

  if (fillnaplot == TRUE) {
    
    colnames <- names(cnmatrix)
    colnames <- colnames[!colnames %in% c("chr", "start", "end", "idx", "width")]
    
    #count proportion of rows that have NA values
    narows <- cnmatrix[, ..colnames]
    narows <- apply(narows, 1, function(x) sum(is.na(x))) / length(colnames)
    #remove the largest region of consecutive NA's, this will the centromere in most chroms
    narowsdf <- cnmatrix[, 1:4] %>% 
      dplyr::mutate(id = 1:dplyr::n()) %>% 
      dplyr::mutate(isna = narows) %>% 
      dplyr::mutate(runid = rleid(isna > 0.9)) %>% 
      dplyr::add_count(runid) %>% 
      dplyr::group_by(chr) %>% 
      dplyr::mutate(maxn = max(n)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(isna & n == maxn & n > 9)
    
    cnmatrix <- cnmatrix %>%
      dplyr::as_tibble() %>%
      tidyr::fill(., colnames, .direction = "downup")
  }
  
  if (fillna == TRUE) {
    colnames <- names(cnmatrix)
    colnames <- colnames[!colnames %in% c("chr", "start", "end", "idx", "width")]
    
    cnmatrix <- cnmatrix %>%
      dplyr::as_tibble() %>%
      tidyr::fill(., colnames, .direction = "updown")
  }

  cnmatrix <- as.data.frame(cnmatrix)

  if (na.rm == TRUE) {
    cnmatrix <- na.omit(cnmatrix)
  }

  rownames(cnmatrix) <- paste(cnmatrix$chr, as.integer(cnmatrix$start), as.integer(cnmatrix$end), sep = "_")
  cnmatrix <- subset(cnmatrix, select = -c(idx))

  if (centromere == TRUE & fillna == TRUE){
    cnmatrix[narowsdf$id,5:ncol(cnmatrix)] <- NA
  }

  return(cnmatrix)
}

#' @export
fixjitter <- function(bps, nextend = 2) {
  x <- as.data.frame(colSums(bps))
  names(x) <- "frequency"
  x$loci <- row.names(x)

  frequency_order <- order(x$frequency, decreasing = TRUE)
  nloci <- dim(bps)[2]

  while (length(frequency_order) > 0) {
    idx <- frequency_order[1]
    minidx <- max(1, idx - nextend)
    maxidx <- min(nloci, idx + nextend)
    temp_mat <- bps[, minidx:maxidx] # extract matrix nextend either side of locus of interest
    bps[, minidx:maxidx] <- 0 # set matrix to 0
    bps[, idx] <- as.numeric(rowSums(temp_mat) > 0)
    frequency_order <- setdiff(frequency_order, minidx:maxidx)
  }

  message(paste0("Original number of loci: ", nloci))
  x <- colSums(bps)
  x <- x[x > 0]
  bps <- bps[, names(x)]
  message(paste0("New number of loci: ", dim(bps)[2]))

  return(as.data.frame(bps))
}

#' @export
createbreakpointmatrix <- function(segs,
                                   transpose = FALSE,
                                   internalonly = TRUE,
                                   use_state = FALSE,
                                   state_remove = 2,
                                   fixjitter = FALSE) {
  options("scipen" = 20)

  if (use_state == FALSE) {
    if (internalonly == TRUE) {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, sep = "_")] %>%
        .[, row_num := .I] %>%
        # add row numbers
        .[, remove_row_num := .I[.N], by = .(cell_id, chr)] %>%
        # find last row in each cell_id - chr group
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
  } else {
    if (internalonly == TRUE) {
      segs_bin <- segs %>%
        as.data.table() %>%
        .[, loci := paste(chr, end - 0.5e6 + 1, end, state, sep = "_")] %>%
        .[, row_num := .I] %>%
        # add row numbers
        .[, remove_row_num := .I[.N], by = .(cell_id, chr)] %>%
        # find last row in each cell_id - chr group
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

  if (transpose == TRUE) {
    segs_mat <- subset(segs_mat, select = -c(loci))
    segs_mat <- t(segs_mat)
    if (fixjitter == TRUE) {
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
getBins <- function(chrom.lengths = hg19_chrlength, binsize = 1e6, chromosomes = NULL) {
  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package \"GenomeInfoDb\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  chrom.lengths <- chrom.lengths[names(chrom.lengths)[names(chrom.lengths) != "M"]]

  chrom.lengths <- chrom.lengths[!is.na(names(chrom.lengths))]
  chrom.lengths <- chrom.lengths[!is.na(chrom.lengths)]
  chroms.in.data <- names(chrom.lengths)
  if (is.null(chromosomes)) {
    chromosomes <- chroms.in.data
  }
  chroms2use <- intersect(chromosomes, chroms.in.data)
  ## Stop if none of the specified chromosomes exist
  if (length(chroms2use) == 0) {
    chrstring <- paste0(chromosomes, collapse = ", ")
    stop("Could not find length information for any of the specified chromosomes: ", chrstring)
  }
  ## Issue warning for non-existent chromosomes
  diff <- setdiff(chromosomes, chroms.in.data)
  if (length(diff) > 0) {
    diffs <- paste0(diff, collapse = ", ")
    warning("Could not find length information for the following chromosomes: ", diffs)
  }

  ### Making fixed-width bins ###

  message("Making fixed-width bins for bin size ", binsize, " ...")
  chrom.lengths.floor <- ceiling(chrom.lengths / binsize) * binsize
  clfloor2use <- chrom.lengths.floor[chroms2use]
  clfloor2use <- clfloor2use[clfloor2use >= binsize]
  if (length(clfloor2use) == 0) {
    stop("All selected chromosomes are smaller than binsize ", binsize)
  }
  bins <- unlist(suppressWarnings(GenomicRanges::tileGenome(clfloor2use, tilewidth = binsize)), use.names = FALSE) #this is a hack for when SVs are at the very end of the chromosome, TODO: change SV location

  suppressWarnings(GenomeInfoDb::seqlevels(bins) <- chroms2use)
  suppressWarnings(GenomeInfoDb::seqlengths(bins) <- chrom.lengths[GenomeInfoDb::seqlevels(bins)])
  skipped.chroms <- setdiff(chromosomes, as.character(unique(GenomeInfoDb::seqnames(bins))))
  suppressWarnings(bins <- GenomeInfoDb::dropSeqlevels(bins, skipped.chroms, pruning.mode = "coarse"))

  if (length(skipped.chroms) > 0) {
    warning("The following chromosomes are smaller than binsize ", binsize, ": ", paste0(skipped.chroms, collapse = ", "))
  }


  bins <- as.data.frame(bins) %>%
    dplyr::select(-strand) %>%
    dplyr::rename(chr = seqnames)

  if (any(bins$width != binsize)) {
    stop("tileGenome failed")
  }

  return(bins)
}

#' @export
widen_haplotypebins <- function(haplotypes, binsize = 5e6) {
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
                       genome_coords = hg19_chrlength) {
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
  if ("state_phase" %in% colnames(CNbins)) {
    widerCNbins <- data.table::as.data.table(CNbinstemp_g[S4Vectors::queryHits(overlaps)]) %>%
      .[, c("strand", "start", "seqnames", "end", "width") := NULL] %>%
      dplyr::bind_cols(., bincoords) %>%
      data.table::setnames(., "seqnames", "chr") %>%
      data.table::setcolorder(., c("chr", "start", "end", "cell_id")) %>%
      .[, .(
        state = as.double(round(median(state, na.rm = TRUE))),
        copy = as.double(median(copy, na.rm = TRUE)),
        A = as.double(floor(median(A))),
        alleleA = sum(alleleA),
        alleleB = sum(alleleB),
        totalcounts = sum(totalcounts)
      ), by = .(chr, start, end, cell_id)] %>%
      .[, BAF := alleleB / totalcounts] %>%
      .[, B := state - A] %>%
      dplyr::select(cell_id, chr, start, end, state, copy, B, A, dplyr::everything()) %>%
      add_states()
  } else {
    widerCNbins <- data.table::as.data.table(CNbinstemp_g[S4Vectors::queryHits(overlaps)]) %>%
      .[, c("strand", "start", "seqnames", "end", "width") := NULL] %>%
      dplyr::bind_cols(., bincoords) %>%
      data.table::setnames(., "seqnames", "chr") %>%
      data.table::setcolorder(., c("chr", "start", "end", "cell_id")) %>%
      .[, .(
        state = round(mean(state, na.rm = TRUE)),
        copy = mean(copy, na.rm = TRUE)
      ), by = .(chr, start, end, cell_id)]
  }

  return(as.data.frame(widerCNbins))
}

#' @export
snv_states <- function(SNV, CNbins) {
  CN <- CNbins %>%
    dplyr::rename(chry = chr, starty = start, endy = end, cell_idy = cell_id) %>%
    as.data.table()

  SNV <- as.data.table(SNV)

  mappedSNVs <- CN[SNV,
    on = .(chry == chr, cell_idy == cell_id, starty < start, endy > start)
  ]
  mappedSNVs <- mappedSNVs %>%
    # .[, end := NULL] %>%
    data.table::setnames(., "chry", "chr") %>%
    data.table::setnames(., "starty", "start") %>%
    data.table::setnames(., "cell_idy", "cell_id") %>%
    data.table::setcolorder(., c("chr", "start", "ref", "alt", "cell_id")) %>%
    .[order(cell_id, chr, start)]

  return(as.data.frame(mappedSNVs))
}

#' @export
map_to_segments <- function(regions, segments){
  segments <- dplyr::mutate(segments, seg_start = start, seg_end = end)
  regions <- dplyr::mutate(regions, pos = start)
  mapped_segs <- as.data.table(regions)[as.data.table(segments),
                   on = .(chr == chr, start > start, start < end)
                   ] %>% 
    .[, start := pos] %>% 
    .[, pos := NULL] %>% 
    .[, start.1 := NULL] %>% 
    na.omit()
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
#' @param full 
#' @param mergesmallarms for very small acrocentric arms just use chromosome
#' @return Character vector, with choromosome arm of given genomic coordinates
#'
#' @export
coord_to_arm <- function(chromosome, 
                         position, 
                         assembly = "hg19", 
                         full = FALSE, 
                         mergesmallarms = FALSE) {
  if (length(chromosome) != length(position)) {
    stop("chromosome and position must have equal length")
  }
  if (!(assembly %in% c("hg38", "hg19", "hg18", "hg17", "hg16"))) {
    stop("Invalid assembly, allowed options are hg38, hg19, hg18, hg17 and hg16")
  }
  if (any(substr(chromosome, 1, 3) != "chr")) {
    chromosome <- paste0("chr", chromosome)
  }
  if (any(!grepl("chr[X-Y]|[0-9]+", chromosome))) {
    stop("Invalid chromosome, must be 1-22, X or Y (or chr1-chr22, chrX or chrY)")
  }
  data("cytoband_map", envir = environment())
  arms <- rep("     ", length(chromosome))

  small <- cytoband_map[[assembly]] %>%
    dplyr::mutate(arm = substr(V4, 1, 1)) %>%
    dplyr::group_by(V1, arm) %>%
    dplyr::summarise(start = min(V2), end = max(V3)) %>%
    dplyr::mutate(width = (end - start) / 1e6) %>%
    dplyr::mutate(small = ifelse(width < 40, TRUE, FALSE)) %>%
    dplyr::group_by(V1) %>%
    dplyr::summarise(small = any(small))

  for (i in unique(chromosome)) {
    map <- cytoband_map[[assembly]][V1 == i]
    arm <- map[(findInterval(position[chromosome == i], map$V3) + 1)]$V4
    if (!full) {
      if (mergesmallarms & small[small$V1 == i, ]$small) {
        arm <- ""
      } else {
        arm <- substr(arm, 1, 1)
      }
    }
    arms[chromosome == i] <- arm
  }

  return(arms)
}

#' @export
create_segments <- function(CNbins, field = "state") {
  newsegs <- CNbins %>%
    data.table::as.data.table() %>%
    .[order(cell_id, chr, start)] %>%
    .[, rlid := data.table::rleid(get(field)), by = cell_id] %>%
    .[, list(
      start = min(start),
      end = max(end)
    ), by = .(cell_id, chr, get(field), rlid)] %>%
    .[order(cell_id, chr, start)] %>%
    dplyr::select(cell_id, chr, start, end, dplyr::everything(), -rlid) %>%
    as.data.frame()
  setnames(newsegs, "get", field)
  return(newsegs)
}

#' @export
create_cntransitions <- function(CNbins, 
                                 field = "state", 
                                 add_orientation = TRUE,
                                 binsize = 0.5e6) {
  
  newsegs <- CNbins %>%
    data.table::as.data.table() %>%
    .[order(cell_id, chr, start)] %>%
    .[, rlid := data.table::rleid(get(field)), by = cell_id] %>%
    .[, list(
      start = min(start),
      end = min(start) + binsize - 1
    ), by = .(cell_id, chr, get(field), rlid)] %>%
    .[order(cell_id, chr, start)] %>%
    dplyr::group_by(cell_id, chr) %>% 
    dplyr::mutate(rown = dplyr::row_number()) %>% 
    dplyr::filter(rown != 1) %>% #filter out first row, count the transition as the bin of the second segment
    dplyr::select(cell_id, chr, start, end, dplyr::everything(), -rlid, -rown) %>%
    as.data.frame()
  setnames(newsegs, "get", field)
  
  if (add_orientation){
    newsegs <- CNbins %>%
      data.table::as.data.table() %>%
      .[order(cell_id, chr, start)] %>%
      .[, rlid := data.table::rleid(get(field)), by = cell_id] %>%
      .[, list(
        start = min(start),
        end = min(start) + binsize - 1
      ), by = .(cell_id, chr, state, rlid)] %>%
      .[order(cell_id, chr, start)] %>%
      dplyr::group_by(cell_id, chr) %>% 
      dplyr::mutate(rown = dplyr::row_number()) %>% 
      dplyr::mutate(state_diff = state - lag(state)) %>% 
      dplyr::filter(rown != 1) %>% #filter out first row, count the transition as the bin of the second segment
      dplyr::mutate(orientation = ifelse(state_diff > 0, "-", "+")) %>% 
      dplyr::select(cell_id, chr, start, end, dplyr::everything(), -rlid, -rown) %>%
      as.data.frame()
  }
  
  return(newsegs)
}

#' @export
orderdf <- function(CNbins) {
  dfchr <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))
  dfchr <- dfchr %>% dplyr::filter(chr %in% unique(CNbins$chr))
  return(CNbins %>%
    as.data.table() %>%
    .[dfchr, on = "chr"] %>%
    .[order(cell_id, idx, start)] %>%
    .[, idx := NULL] %>%
    as.data.frame())
}

densmode <- function(x) {
  dens <- density(x)
  dens$x[which.max(dens$y)]
}

#' @export
qc_summary <- function(cn) {
  if (is.hscn(cn) | is.ascn(cn)) {
    cn <- cn$data
  } else {
    cn <- cn
  }

  distance_df <- cn %>%
    dplyr::filter(state_AS_phased != "0|0") %>%
    dplyr::group_by(state_AS_phased, B, A) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frac = n / dplyr::n()) %>%
    dplyr::filter(n > 10) %>%
    dplyr::group_by(state_AS_phased, B, A, n, frac) %>%
    dplyr::summarise(
      medianBAF = median(BAF, na.rm = TRUE),
      meanBAF = mean(BAF, na.rm = TRUE),
      # modeBAF = densmode(BAF),
      high95 = quantile(BAF, 0.975, na.rm = TRUE),
      low95 = quantile(BAF, 0.025, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(expBAF = B / (B + A)) %>%
    dplyr::mutate(distance = sqrt((expBAF - medianBAF)^2)) %>%
    as.data.frame()

  summary_distance <- weighted.mean(distance_df$distance, distance_df$frac)

  message(paste0("Average distance from median to expected BAF = ", round(summary_distance, 4)))

  return(list(distance = distance_df, summary = summary_distance))
}

#' @export
qc_per_cell <- function(cn){
  
  if (is.hscn(cn) | is.ascn(cn)) {
    cn <- cn$data
  } else {
    cn <- cn
  }
  
  qc <- cn %>%
    dplyr::filter(state_AS_phased != "0|0") %>%
    dplyr::group_by(state_AS_phased, B, A, cell_id) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(frac = n / dplyr::n()) %>%
    dplyr::filter(n > 10) %>%
    dplyr::group_by(state_AS_phased, B, A, n, frac, cell_id) %>%
    dplyr::summarise(
      medianBAF = median(BAF, na.rm = TRUE),
      meanBAF = mean(BAF, na.rm = TRUE),
      # modeBAF = densmode(BAF),
      high95 = quantile(BAF, 0.975, na.rm = TRUE),
      low95 = quantile(BAF, 0.025, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(expBAF = B / (B + A)) %>%
    dplyr::mutate(distance = sqrt((expBAF - medianBAF)^2)) %>%
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(average_distance = weighted.mean(distance, frac)) %>% 
    as.data.frame(.)
  
  pl <- cn %>% as.data.table() %>% .[!chr %in% c("X", "Y"), 
                                     list(ploidy = Mode(state), totalhapcounts = sum(totalcounts)), 
                                     by = "cell_id"]
  
  nsegs <- create_segments(cn) %>% 
    dplyr::filter(!chr %in% c("X", "Y")) %>% 
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(nsegments = dplyr::n() - 22)
  
  qc <- dplyr::left_join(qc, pl, by = "cell_id") %>% 
    dplyr::left_join(., nsegs, by = "cell_id") %>% 
    as.data.frame()
  return(qc)
}

#' @export
per_arm_baf_mat <- function(haps,
                            mergelowcounts = TRUE,
                            mincounts = 10,
                            arms = NULL) {
  if (mergelowcounts & is.null(arms)) {
    message(paste0("Only using arms with at least an average ", mincounts, " counts"))
    baf <- haps %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        BAF = alleleB / (alleleA + alleleB),
        total = alleleB + alleleA
      )

    mergearms <- baf %>%
      dplyr::group_by(chr, chrarm, arm) %>%
      dplyr::summarise(median = median(total)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(merge = ifelse(median < mincounts, "merge", "no")) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(merge = ifelse(any(merge == "merge"), "merge", "no"))

    baf <- haps %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::left_join(mergearms) %>%
      dplyr::mutate(arm = ifelse(merge == "merge", "", arm)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        BAF = alleleB / (alleleA + alleleB),
        total = alleleB + alleleA
      )
  } else if (mergelowcounts == FALSE & is.null(arms)) {
    message("Using all chromosome arms")
    baf <- haps %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        BAF = alleleB / (alleleA + alleleB),
        total = alleleB + alleleA
      )
  } else if (!is.null(arms)) {
    message(paste0("Only using specific chromosome arms: "), paste0(arms, collapse = ", "))
    baf <- haps %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(arm = ifelse(chrarm %in% arms, arm, "")) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        BAF = alleleB / (alleleA + alleleB),
        total = alleleB + alleleA
      )
  }

  ord <- dplyr::distinct(baf, chr, arm, chrarm) %>%
    as.data.table() %>%
    .[gtools::mixedorder(chrarm)] %>%
    .[, idx := 1:.N] %>%
    dplyr::as_tibble()

  baf <- dplyr::left_join(baf, ord)

  baf_mat <- baf %>%
    dplyr::select(cell_id, chrarm, BAF) %>%
    tidyr::pivot_wider(names_from = "chrarm", values_from = "BAF") %>%
    as.data.frame()

  row.names(baf_mat) <- baf_mat$cell_id
  baf_mat <- subset(baf_mat, select = -c(cell_id))

  idx <- dplyr::distinct(baf, chr, arm, chrarm, idx) %>% dplyr::arrange(idx)
  baf_mat <- baf_mat[, idx$chrarm]

  counts <- baf %>%
    dplyr::group_by(chr, chrarm, arm) %>%
    dplyr::summarise(median = median(total)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(median < mincounts)

  warning(paste0("The following chromosomes have on average ", mincounts, " or fewer counts: ", paste(counts$chrarm, collapse = ", ")))

  return(list(bafperchr = baf, bafperchrmat = baf_mat))
}

mixedrank = function(x) order(gtools::mixedorder(x))

#' @export
filter_segments <- function(segs, binwidth = 5e6){
  data("hg19chrom_coordinates", envir = environment())
  chrin <- unique(segs$chr)
  
  w <- 0
  while (w < 5e6){
    segs <- segs %>% 
      dplyr::mutate(w = end-start) %>% 
      dplyr::filter(w > binwidth) %>% 
      dplyr::group_by(chr, cell_id) %>% 
      dplyr::mutate(end2 = dplyr::lead(start) - 1) %>% 
      dplyr::mutate(end = ifelse(is.na(end2), end, end2)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(chr, start, end, w) %>% 
      dplyr::arrange(mixedrank(chr), start)
    w <- min(segs$w)
  }
  
  segs <- dplyr::select(segs, -w)
  
  #make sure whole chromosome is covered
  segs <- dplyr::left_join(segs, hg19chrom_coordinates %>% dplyr::filter(arm == "") %>% dplyr::mutate(start = start + 1), by = "chr") %>% 
    dplyr::group_by(chr) %>%  
    dplyr::mutate(start = ifelse(dplyr::row_number() == 1, start.y, start.x), 
                  end = ifelse(dplyr::row_number() == dplyr::n(), end.y, end.x)) %>% 
    dplyr::select(chr, start, end)
  
  w <- 0
  while (w < 5e6){
    segs <- segs %>% 
      dplyr::mutate(w = end-start) %>% 
      dplyr::filter(w > binwidth) %>% 
      dplyr::group_by(chr) %>% 
      dplyr::mutate(end2 = dplyr::lead(start) - 1) %>% 
      dplyr::mutate(end = ifelse(is.na(end2), end, end2)) %>% 
      dplyr::ungroup() %>% 
      dplyr::select(chr, start, end, w) %>% 
      dplyr::arrange(mixedrank(chr), start)
    w <- min(segs$w)
  }
  
  segs <- dplyr::select(segs, -w)
  
  if (all(chrin %in% unique(segs$chr)) == FALSE){
    warning("Not all chromosomes present in output adding whole chromosomes")
    chroms <- chrin[which(!chrin %in% unique(segs$chr))]
    to_add <- hg19chrom_coordinates %>% 
      dplyr::filter(chr %in% chroms & arm == "") %>% 
      dplyr::mutate(start = start + 1) %>% 
      dplyr::select(chr, start, end)
    segs <- dplyr::bind_rows(segs, to_add) %>% 
      dplyr::arrange(mixedrank(chr), start)
  }
    
  return(segs %>% as.data.frame())
}

#' @export
per_segment_baf_mat <- function(haps,
                            segments = NULL,
                            mergelowcounts = TRUE,
                            mincounts = 10,
                            arms = NULL) {
  
  options(scipen=999)
  
  if (is.null(segments)){
    message("No segments provided, using chromosome arms as segments")
    segments <- hg19chrom_coordinates %>% 
      dplyr::filter(arm != "") %>% 
      dplyr::select(-arm, -chrarm) %>% 
      dplyr::mutate(start = start + 1)
  }
  
  mapped_segs <- map_to_segments(haps, segments)
  
  mapped_segs$segid <- paste(mapped_segs$chr, mapped_segs$seg_start, mapped_segs$seg_end, sep = "_")
  
  baf <- mapped_segs %>%
    dplyr::group_by(chr, segid, cell_id) %>%
    dplyr::summarise(
      alleleA = sum(alleleA, na.rm = TRUE),
      alleleB = sum(alleleB, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      BAF = alleleB / (alleleA + alleleB),
      total = alleleB + alleleA
    )
  
  ord <- dplyr::distinct(baf, chr, segid) %>%
    as.data.table() %>%
    .[gtools::mixedorder(chr)] %>%
    .[, idx := 1:.N] %>%
    dplyr::as_tibble()
  
  baf <- dplyr::left_join(baf, ord)
  
  baf_mat <- baf %>%
    dplyr::select(cell_id, segid, BAF) %>%
    tidyr::pivot_wider(names_from = "segid", values_from = "BAF") %>%
    as.data.frame()
  
  row.names(baf_mat) <- baf_mat$cell_id
  baf_mat <- subset(baf_mat, select = -c(cell_id))
  
  idx <- dplyr::distinct(baf, chr, segid, idx) %>% dplyr::arrange(idx)
  baf_mat <- baf_mat[, idx$segid]
  
  counts <- baf %>%
    dplyr::group_by(chr, segid) %>%
    dplyr::summarise(median = median(total)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(median < mincounts)
  
  warning(paste0("The following segments have on average ", mincounts, " or fewer counts: ", paste(counts$segid, collapse = ", ")))
  
  return(list(bafpersegment = baf, bafpersegmentmat = baf_mat))
}

#' @export
per_chrarm_cn <- function(hscn, arms = NULL) {
  data("hg19chrom_coordinates", envir = environment())

  if ("state_phase" %in% colnames(hscn)){
    if (is.null(arms)) {
      hscn_arm <- hscn %>%
        dplyr::filter(chr != "Y") %>%
        dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        as.data.table() %>%
        filter(!chrarm %in% c("13p", "14p", "15p", "21p", "22p", "Yp")) %>% #remove acrocentric chromosomes
        .[, list(
          state = as.double(round(median(state, na.rm = TRUE))),
          copy = as.double(median(copy, na.rm = TRUE)),
          A = as.double(floor(median(A, na.rm = TRUE))),
          alleleA = sum(alleleA, na.rm = TRUE),
          alleleB = sum(alleleB, na.rm = TRUE),
          totalcounts = sum(totalcounts, na.rm = TRUE),
          state_sd = sd(state, na.rm = TRUE),
          proportion = sum(state_AS_phased == Mode(state_AS_phased)) / .N
        ), by = c("chr", "arm", "chrarm", "cell_id")] %>%
        .[, BAF := alleleB / totalcounts] %>%
        .[, B := state - A]
    } else {
      hscn_arm <- hscn %>%
        dplyr::filter(chr != "Y") %>%
        dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        dplyr::mutate(arm = ifelse(chrarm %in% arms, arm, "")) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        filter(!chrarm %in% c("13p", "14p", "15p", "21p", "22p", "Yp")) %>%
        as.data.table() %>%
        .[, list(
          state = as.double(round(median(state, na.rm = TRUE))),
          copy = as.double(median(copy, na.rm = TRUE)),
          A = as.double(floor(median(A, na.rm = TRUE))),
          alleleA = sum(alleleA), na.rm = TRUE,
          alleleB = sum(alleleB, na.rm = TRUE),
          totalcounts = sum(totalcounts, na.rm = TRUE),
          state_sd = sd(state, na.rm = TRUE),
          proportion = sum(state_AS_phased == Mode(state_AS_phased)) / .N
        ), by = c("chr", "arm", "chrarm", "cell_id")] %>%
        .[, BAF := alleleB / totalcounts] %>%
        .[, B := state - A]
    }
  
    hscn_arm <- hscn_arm %>%
      add_states() %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      dplyr::mutate(state = ifelse(state > 11, 11, state))
    
  } else {
    if (is.null(arms)) {
      hscn_arm <- hscn %>%
        dplyr::filter(chr != "Y") %>%
        dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        as.data.table() %>%
        .[, list(
          state = as.double(round(median(state, na.rm = TRUE))),
          copy = as.double(median(copy, na.rm = TRUE)),
          state_sd = sd(state, na.rm = TRUE),
          proportion = sum(state== Mode(state)) / .N
        ), by = c("chr", "arm", "chrarm", "cell_id")]
    } else {
      hscn_arm <- hscn %>%
        dplyr::filter(chr != "Y") %>%
        dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        dplyr::mutate(arm = ifelse(chrarm %in% arms, arm, "")) %>%
        dplyr::mutate(chrarm = paste0(chr, arm)) %>%
        as.data.table() %>%
        .[, list(
          state = as.double(round(median(state, na.rm = TRUE))),
          copy = as.double(median(copy, na.rm = TRUE)),
          state_sd = sd(state, na.rm = TRUE),
          proportion = sum(state == Mode(state)) / .N
        ), by = c("chr", "arm", "chrarm", "cell_id")] 
    }
    
    hscn_arm <- hscn_arm %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      dplyr::mutate(state = ifelse(state > 11, 11, state))
  }

  return(hscn_arm)
}

#' @export
per_chr_cn <- function(hscn, arms = NULL) {
  data("hg19chrom_coordinates", envir = environment())
  
  if ("state_phase" %in% colnames(hscn)){

    hscn_chr <- hscn %>%
      dplyr::filter(chr != "Y") %>%
      as.data.table() %>%
      .[, list(
        state = as.double(round(median(state, na.rm = TRUE))),
        copy = as.double(median(copy, na.rm = TRUE)),
        A = as.double(floor(median(A, na.rm = TRUE))),
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE),
        totalcounts = sum(totalcounts, na.rm = TRUE),
        state_sd = sd(state, na.rm = TRUE),
        proportion = sum(state_AS_phased == Mode(state_AS_phased)) / .N
      ), by = c("chr", "cell_id")] %>%
      .[, BAF := alleleB / totalcounts] %>%
      .[, B := state - A] %>%
      add_states() %>%
      dplyr::left_join(hg19chrom_coordinates %>% dplyr::filter(arm == "")) %>%
      dplyr::mutate(state = ifelse(state > 11, 11, state))
  } else {
    hscn_chr <- hscn %>%
      dplyr::filter(chr != "Y") %>%
      as.data.table() %>%
      .[, list(
        state = as.double(round(median(state, na.rm = TRUE))),
        copy = as.double(median(copy, na.rm = TRUE)),
        state_sd = sd(state, na.rm = TRUE),
        proportion = sum(state == Mode(state)) / .N
      ), by = c("chr", "cell_id")] %>%
      dplyr::left_join(hg19chrom_coordinates %>% dplyr::filter(arm == "")) %>%
      dplyr::mutate(state = ifelse(state > 11, 11, state))
  }

  return(hscn_chr)
}

#' @export
add_states <- function(df) {
  df <- df %>%
    .[, state_AS_phased := paste0(A, "|", B)] %>%
    .[, state_AS := paste0(pmax(state - B, B), "|", pmin(state - B, B))] %>%
    .[, state_min := pmin(A, B)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
      1 * ((B < A)) +
      2 * ((B > A))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-Hom", "B-Hom")[1 +
      1 * ((B < A) & (B != 0)) +
      2 * ((B > A) & (A != 0)) +
      3 * ((B < A) & (B == 0)) +
      4 * ((B > A) & (A == 0))]] %>%
    .[, state_BAF := round((B / state) / 0.1) * 0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]
  return(df)
}

#' @export
createBAFassay <- function(seur, rna_ascn, ref = "hg19") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package \"Seurat\" is needed for this function. Please install it.",
      call. = FALSE
    )
  }
  
  data("gene_locations", envir = environment())

  message("Add BAF to Seurat object")
  x <- tidyr::pivot_wider(rna_ascn$hscn %>%
    dplyr::select(cell_id, BAF, segid) %>%
    dplyr::mutate(segid = paste0("BAF-", stringr::str_replace_all(segid, "_", "-"))),
  names_from = "segid",
  values_from = c("BAF")
  ) %>%
    as.data.frame()

  row.names(x) <- x$cell_id
  x <- as.matrix(subset(x, select = -cell_id))
  cM <- colMeans(x, na.rm = TRUE)
  indx <- which(is.na(x), arr.ind = TRUE)
  x[indx] <- cM[indx[, 2]]
  baf_assay <- Seurat::CreateAssayObject(data = t(x))
  seur[["segBAF"]] <- Seurat::CreateAssayObject(data = t(x))

  message("Add allele specific state to Seurat Object")
  x <- tidyr::pivot_wider(rna_ascn$hscn %>%
    dplyr::select(cell_id, state_phase, segid),
  names_from = "segid",
  values_from = c("state_phase")
  ) %>%
    as.data.frame()
  row.names(x) <- x$cell_id
  x <- as.matrix(subset(x, select = -cell_id))
  seur[["ASDP"]] <- Seurat::CreateAssayObject(data = t(x))

  message("Add clone id to metadata")
  clonesdf <- dplyr::distinct(rna_ascn$clusters, cell_id, clone_id) %>%
    as.data.frame() %>%
    dplyr::mutate(clone_id = ifelse(clone_id == "" | is.na(clone_id), NA, clone_id)) %>%
    tidyr::fill(clone_id, .direction = "downup")
  clonesvec <- clonesdf$clone_id
  names(clonesvec) <- clonesdf$cell_id
  seur <- AddMetaData(
    object = seur,
    metadata = clonesvec,
    col.name = "DP_cloneid"
  )
  
  if (!is.null(rna_ascn$haplotypecounts)){
    message("Add gene allelic imbalance")
    snps_to_genes <- map_to_segments(rna_ascn$haplotypecounts, 
                                     gene_locations[[ref]] %>% dplyr::filter(!stringr::str_detect(chr, "_")))
    gene_baf <- snps_to_genes %>% 
      .[, list(BAF = mean(BAF)), by = c("cell_id", "ensembl_gene_symbol")] %>% 
      .[, cell_id := as.factor(cell_id)] %>% 
      .[, ensembl_gene_symbol := as.factor(paste0("BAF-",ensembl_gene_symbol))]
    x <- tidyr::pivot_wider(gene_baf %>%
                              dplyr::select(cell_id, ensembl_gene_symbol, BAF),
                            names_from = "ensembl_gene_symbol",
                            values_from = c("BAF")
    ) %>%
      as.data.frame()
    row.names(x) <- x$cell_id
    x <- as.matrix(subset(x, select = -cell_id))
    seur[["gBAF"]] <- Seurat::CreateAssayObject(data = t(x))
    meta <- gene_locations[[ref]] %>% 
      dplyr::filter(!stringr::str_detect(chr, "_")) %>% 
      dplyr::group_by(ensembl_gene_symbol) %>% 
      dplyr::filter(dplyr::row_number() == 1) %>%  #hack, take first chr if there are > 1
      as.data.frame() %>% 
      dplyr::mutate(arm = paste0(chr, coord_to_arm(chr, start, assembly = ref)))
    row.names(meta) <- meta$ensembl_gene_symbol
    meta <- meta[stringr::str_remove(rownames(seur[["gBAF"]]), "BAF-"), ]
    meta <- merge(meta, snps_to_genes[, list(totalsnpcounts = sum(totalcounts)), by = "ensembl_gene_symbol"])
    row.names(meta) <- meta$ensembl_gene_symbol
    seur[["gBAF"]]@meta.features <- meta
    
    gene_counts <- snps_to_genes %>% 
      .[, list(totalcounts = sum(totalcounts)), by = c("cell_id", "ensembl_gene_symbol")] %>% 
      .[, cell_id := as.factor(cell_id)] %>% 
      .[, ensembl_gene_symbol := as.factor(paste0("counts-",ensembl_gene_symbol))]
    x <- tidyr::pivot_wider(gene_counts %>%
                              dplyr::select(cell_id, ensembl_gene_symbol, totalcounts),
                            names_from = "ensembl_gene_symbol",
                            values_from = c("totalcounts")
    ) %>%
      as.data.frame()
    row.names(x) <- x$cell_id
    x <- as.matrix(subset(x, select = -cell_id))
    seur[["gCounts"]] <- Seurat::CreateAssayObject(data = t(x))
  }
  
  # x <- with(gene_baf, Matrix::sparseMatrix(i=as.numeric(ensembl_gene_symbol),
  #                                          j=as.numeric(cell_id), 
  #                                          x=as.numeric(BAF),
  #                                          dimnames=list(levels(ensembl_gene_symbol), levels(cell_id))))
  #seur[["gBAF"]] <- Seurat::CreateAssayObject(data = x)

  return(seur)
}

#' Add gene genomic location to seurat object
#' 
#' @param obj Seurat object
#' @param ref Human genome reference, default `hg19`
#' 
#' export
add_gene_locations_to_seurat <- function(obj, ref = "hg19"){
  meta <- gene_locations[[ref]] %>% 
    dplyr::filter(!stringr::str_detect(chr, "_")) %>% 
    dplyr::group_by(ensembl_gene_symbol) %>% 
    dplyr::filter(dplyr::row_number() == 1) %>%  #hack, take first chr if there are > 1
    as.data.frame() %>% 
    dplyr::mutate(arm = paste0(chr, coord_to_arm(chr, start, assembly = ref)))
  row.names(meta) <- meta$ensembl_gene_symbol
  meta <- meta[rownames(obj[["RNA"]]), ]
  row.names(meta) <- row.names(obj[["RNA"]]@meta.features)
  obj[["RNA"]]@meta.features <- cbind(obj[["RNA"]]@meta.features, meta)
  return(obj)
}

#' @export
consensuscopynumber <- function(hscn, cl = NULL) {
  if (!is.null(cl)) {
    hscn <- hscn %>% dplyr::filter(cell_id %in% cl$cell_id)
    hscn <- dplyr::left_join(hscn, cl,  by = "cell_id")
  } else {
    hscn$clone_id <- "Merged Cells"
  }

  if ("state_phase" %in% colnames(hscn)) {
    cn <- hscn %>%
      as.data.table(.) %>%
      .[, .(
        state = as.double(round(median(state, na.rm = TRUE))),
        copy = as.double(median(copy, na.rm = TRUE)),
        A = as.double(floor(median(A, na.rm = TRUE))),
        alleleA = sum(alleleA, na.rm = TRUE),
        alleleB = sum(alleleB, na.rm = TRUE),
        totalcounts = sum(totalcounts, na.rm = TRUE),
        BAF = median(BAF, na.rm = TRUE)
      ), by = .(chr, start, end, clone_id)] %>%
      .[, A := ifelse(A > state, state, A)] %>% 
      .[, B := state - A] %>%
      add_states() %>%
      dplyr::select(clone_id, chr, start, end, state, copy, B, A, dplyr::everything()) %>%
      dplyr::rename(cell_id = clone_id) %>%
      as.data.frame(.)
  } else {
    cn <- hscn %>%
      as.data.table(.) %>%
      .[, .(
        state = round(median(state, na.rm = TRUE)),
        copy = median(copy, na.rm = TRUE)
      ), by = .(chr, start, end, clone_id)] %>%
      dplyr::rename(cell_id = clone_id) %>%
      as.data.frame(.)
  }

  return(cn)
}

#' @export
singletons <- function(cn, field) {
  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }

  x <- CNbins %>%
    dplyr::mutate(changepoint_forward = dplyr::lag({{ field }}) != {{ field }}) %>%
    dplyr::mutate(changepoint_reverse = dplyr::lead({{ field }}) != {{ field }}) %>%
    dplyr::mutate(singleton = changepoint_forward & changepoint_reverse) %>%
    dplyr::select(-changepoint_forward, -changepoint_reverse)

  x_s <- sum(x$singleton)
  message(paste0("Number of singletons: ", x_s))

  return(x)
}

#' @export
BAFdistance <- function(cn) {
  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }

  celldist <- cn %>%
    as.data.table() %>%
    .[, dist := sqrt((BAF - (B / state))^2)] %>%
    .[, dist := BAF - (B / state)] %>%
    .[, list(BAF_distance = mean(dist, na.rm = T)), by = c("cell_id", "chr")]

  return(celldist)
}

#' @export
filtercn <- function(cn, 
                     average_distance_filt = 0.1, 
                     ploidy_filt = NULL, 
                     totalhapcounts_filt = 0, 
                     nsegments_filt = 0){
  qcnew <- dplyr::filter(cn$qc_per_cell, average_distance <= average_distance_filt,
                             totalhapcounts >= totalhapcounts_filt,
                              nsegments >= nsegments_filt)
  if (!is.null(ploidy_filt)){
    qcnew <- dplyr::filter(qcnew, ploidy == ploidy_filt)
  }
  
  cells <- qcnew$cell_id
  cn$data <- cn$data %>% 
    dplyr::filter(cell_id %in% cells)
  cn$qc_summary <- qc_summary(cn)
  cn$qc_per_cell <- qcnew
  
  return(cn)
}

#' @export
filterbycells <- function(cn, cells){
  cn$data <- cn$data %>% 
    dplyr::filter(cell_id %in% cells)
  cn$qc_per_cell <- cn$qc_per_cell %>% 
    dplyr::filter(cell_id %in% cells)
  cn$qc_summary <- qc_summary(cn)
  
  return(cn)
}

#' @export
format_tree_labels <- function(tree, removeloci = TRUE, internal_node_string = 'locus|internal'){
  if (removeloci){
    tip.loci <- grep(internal_node_string, tree$tip.label, value = T)
    while (length(tip.loci) > 0) {
      tree <- ape::drop.tip(tree, tip.loci, trim.internal = FALSE, collapse.singles = FALSE)
      tip.loci <- grep(internal_node_string, tree$tip.label, value = T)
    }
  }
  tree$tip.label <- str_remove(tree$tip.label, "cell_")
  return(tree)
}
