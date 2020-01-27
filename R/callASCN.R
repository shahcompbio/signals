#' @export
combineBAFCN <- function(haplotypes, CNbins, binsize = 0.5e6, filtern = 0){

  message("Finding overlapping cell IDs between CN data and haplotype data...")
  cellidoverlap <- intersect(CNbins$cell_id, haplotypes$cell_id)

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
  haplotypes <- haplotypes %>%
    format_haplotypes(.)
  haplotypes <- data.table::as.data.table(haplotypes)

  if (binsize > 0.5e6){
    message("Widening CN bins")
    CNbins <- CNbins %>%
      widen_bins(., binsize = binsize) %>%
      data.table::as.data.table(.)
    message("Widening haplotype bins")
    haplotypes <- haplotypes %>%
      widen_haplotypebins(., binsize = binsize) %>%
      data.table::as.data.table(.)
  }

  message("Joining bins and haplotypes...")
  #CNbins <- data.table::merge.data.table(CNbins, haplotypes)
  CNbins <- CNbins[haplotypes, on = c("chr", "start", "end", "cell_id"), nomatch=0]

  message("Calculate BAF per bin...")
  CNBAF <- CNbins %>%
    data.table::as.data.table() %>%
    .[totalcounts > filtern] %>%
    .[, lapply(.SD, sum), by = .(chr, start, end, cell_id, state, copy, reads), .SDcols = c("alleleA", "alleleB", "totalcounts")] %>%
    .[, BAF := alleleB / totalcounts]

  return(as.data.frame(CNBAF))
}

#' @export
callAlleleSpecificCN <- function(CNBAF, maxCN = 11){
  maj <- seq(0, maxCN, 1)
  min <- seq(0, maxCN, 1)
  allASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS_phased = paste0(state-min, "|", min)) %>%
    dplyr::mutate(Maj = state-min, Min = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(Maj >= 0, Min >= 0) %>%
    data.table::as.data.table()
  allASstates$cBAF[is.nan(allASstates$cBAF)] <- 0.0

  CNBAF <- data.table::as.data.table(CNBAF)
  CNBAF$copy[is.nan(CNBAF$copy)] <- NA
  setnafill(CNBAF, type = "locf", cols = c("copy"))
  setnafill(CNBAF, type = "nocb", cols = c("copy"))

  alleleCN <- data.table::merge.data.table(CNBAF, allASstates, by = .EACHI, allow.cartesian=TRUE) %>%
    .[, dist := sqrt((copy - state)^2 + (BAF - cBAF)^2)]

  data.table::setcolorder(alleleCN, c("chr", "start", "end", "cell_id", "state"))

  alleleCN <- alleleCN[alleleCN[, .I[which.min(dist)], by = .(chr, start, cell_id)]$V1]

  alleleCN <- alleleCN %>%
    .[, c("cBAF", "dist") := NULL] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                                 1 * ((Min < Maj) & (Min != 0)) +
                                                                                 2 * ((Min > Maj) & (Maj != 0)) +
                                                                                 3 * ((Min < Maj) & (Min == 0)) +
                                                                                 4 * ((Min > Maj) & (Maj == 0))]
      ] %>%
    .[, c("Maj", "Min") := NULL] %>%
    .[order(cell_id, chr, start)]

  return(as.data.frame(alleleCN))
}

#' @export
alleleHMM <- function(n,
                      x,
                      binstates,
                      minor_cn,
                      loherror = 0.01,
                      selftransitionprob = 0.999,
                      eps = 1e-12){

  minor_cn_mat <- t(replicate(length(binstates), minor_cn))
  total_cn_mat <- replicate(length(minor_cn), binstates)

  p <- t(vapply(binstates, function(x) minor_cn / x, FUN.VALUE = numeric(length(minor_cn))))
  p[, 1] <- loherror
  p[minor_cn_mat == total_cn_mat] <- 1 - loherror

  l <- suppressWarnings(dbinom(x, n, p, log = T))
  if (eps > 0.0){
    l <- matrix(vapply(l, function(x) matrixStats::logSumExp(c(x, log(eps))), FUN.VALUE = numeric(1)), dim(l)[1], dim(l)[2])
  }
  l[is.na(l)] <- log(0.0)
  l[minor_cn_mat > total_cn_mat] <- log(0.0)

  tProbs <- matrix((1 - selftransitionprob) / (length(minor_cn) -1), length(minor_cn), length(minor_cn))
  diag(tProbs) <- selftransitionprob
  colnames(tProbs) <- paste0(minor_cn)
  row.names(tProbs) <- paste0(minor_cn)

  states <- paste0(minor_cn)

  res <- myviterbi(l, log(tProbs), observations = 1:length(binstates))

  return(list(minorcn = res, l = l))
}

#' @export
assignalleleHMM <- function(CNBAF,
                            minor_cn,
                            eps = 1e-12,
                            loherror = 0.01,
                            selftransitionprob = 0.999,
                            pb = NULL){

  if (!is.null(pb)){
    pb$tick()$print()
  }

  hmmresults <- alleleHMM(n = CNBAF$totalcounts,
                          x = CNBAF$alleleB,
                          CNBAF$state,
                          minor_cn,
                          loherror = loherror,
                          eps = eps,
                          selftransitionprob = selftransitionprob)

  CNBAF$state_min <- as.numeric(hmmresults$minorcn)

  CNBAF <- data.table::as.data.table(CNBAF)

  CNBAF <- CNBAF %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
                                           1 * ((Min < Maj)) +
                                           2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                               1 * ((Min < Maj) & (Min != 0)) +
                                                                               2 * ((Min > Maj) & (Maj != 0)) +
                                                                               3 * ((Min < Maj) & (Min == 0)) +
                                                                               4 * ((Min > Maj) & (Maj == 0))]
    ] %>%
  #.[, c("Maj", "Min") := NULL] %>%
  .[order(cell_id, chr, start)]

  return(as.data.frame(CNBAF))
}

#' @export
callalleleHMMcell <- function(CNBAF,
                              minor_cn,
                              eps = 1e-9,
                              loherror = 0.01,
                              selftransitionprob = 0.999){

  hmmresults <- alleleHMM(n = CNBAF$totalcounts,
                          x = CNBAF$alleleB,
                          CNBAF$state,
                          minor_cn,
                          loherror = loherror,
                          eps = eps,
                          selftransitionprob = selftransitionprob)

  CNBAF$state_min <- as.numeric(hmmresults$minorcn)
  CNBAF$allele <- ifelse()

  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
                                         1 * ((Min < Maj)) +
                                         2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                                 1 * ((Min < Maj) & (Min != 0)) +
                                                                                 2 * ((Min > Maj) & (Maj != 0)) +
                                                                                 3 * ((Min < Maj) & (Min == 0)) +
                                                                                 4 * ((Min > Maj) & (Maj == 0))]
      ] %>%
    #.[, c("Maj", "Min") := NULL] %>%
    .[order(cell_id, chr, start)]

  return(list(alleleCN = CNBAF, posterior_prob = hmmresults$posterior_prob, l = hmmresults$l))
}

#' @export
callAlleleSpecificCNHMM <- function(CNBAF,
                                    eps = 1e-12,
                                    loherror = 0.01,
                                    maxCN = 10,
                                    selftransitionprob = 0.999,
                                    progressbar = TRUE,
                                    ncores = 1){

  minor_cn <- seq(0, maxCN, 1)

  if (progressbar == TRUE){
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else{
    pb <- NULL
  }

  #Make sure dataframe is in chromosome position order
  CNBAF <- CNBAF %>%
    data.table::as.data.table() %>%
    .[order(cell_id, chr, start)]

  if (ncores > 1){
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
                                                    function(cell) assignalleleHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                                   eps = eps, loherror = loherror,
                                                                                   selftransitionprob = selftransitionprob,
                                                                                   pb = pb), mc.cores = ncores)) %>%
      .[order(cell_id, chr, start)]
  } else{
    alleleCN <- data.table::rbindlist(lapply(unique(CNBAF$cell_id),
                                        function(cell) assignalleleHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                       eps = eps, loherror = loherror,
                                                                       selftransitionprob = selftransitionprob,
                                                                       pb = pb))) %>%
      .[order(cell_id, chr, start)]
  }

  return(as.data.frame(alleleCN))
}

#' @export
countsegments <- function(ascn, armlevel = TRUE, binfilter = 10){
  segments <- create_segments(ascn, phase)
  segments$chrarm <- paste0(segments$chr, coord_to_arm(segments$chr, segments$start))
  if (armlevel == TRUE){
    segments_counts <- segments %>%
      dplyr::filter(nbin >= binfilter) %>%
      dplyr::group_by(cell_id, chrarm, phase) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::pivot_wider(values_from = n, names_from = phase) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(entB = A * log(A) + B * log(B)) %>%
      dplyr::arrange(desc(entB)) %>%
      dplyr::group_by(chrarm) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup()
  } else {
    segments_counts <- segments %>%
      dplyr::filter(nbin >= binfilter) %>%
      dplyr::group_by(cell_id, chr, phase) %>%
      dplyr::summarise(n = n()) %>%
      tidyr::pivot_wider(values_from = n, names_from = phase) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(entB = A * log(A) + B * log(B)) %>%
      dplyr::arrange(desc(entB)) %>%
      dplyr::group_by(chr) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup()
  }
  return(segments_counts)
}

#' @export
countchangesegments <- function(ascn, armlevel = TRUE, binfilter = 10){
  segments <- create_segments(ascn, phase)
  segments_state <- create_segments(ascn, state)
  percell <- segments_state %>%
    dplyr::group_by(chr, cell_id) %>%
    dplyr::summarise(nstates = n()) %>%
    dplyr::ungroup()
  segments$chrarm <- paste0(segments$chr, coord_to_arm(segments$chr, segments$start))
  if (armlevel == TRUE){
    segments_counts <- segments %>%
      dplyr::filter(nbin >= binfilter, phase != "Balanced") %>%
      dplyr::group_by(cell_id, chrarm) %>%
      dplyr::mutate(change = phase != lag(phase)) %>%
      dplyr::summarise(n = sum(change, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::group_by(chrarm) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup()
  } else {
    segments_counts <- segments %>%
      dplyr::filter(nbin >= binfilter) %>%
      dplyr::group_by(cell_id, chr) %>%
      dplyr::mutate(change = (phase != lag(phase)) & (phase != "Balanced")) %>%
      dplyr::summarise(n = sum(change, na.rm = T)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(percell) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(filtn = max(5, min(nstates))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(nstates <= filtn) %>%
      dplyr::arrange(desc(n)) %>%
      dplyr::group_by(chr) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::ungroup()
  }
  return(segments_counts)
}


find_switches <- function(ascn, armlevel = FALSE, binfilter = 10){
  bins <- dplyr::distinct(ascn, chr, start, end)
  cell_segments <- countchangesegments(ascn, armlevel = armlevel, binfilter = binfilter)
  bin_switch <- data.frame()
  for (i in 1:nrow(cell_segments)){
    print(i)
    print(cell_segments$chr[i])
    print(ascn %>%
            filter(chr == cell_segments$chr[i]) %>%
            distinct(chr, start, end) %>%
            dim())
    if (cell_segments$n[i] < 2){
      print("filtered")
      ascn_temp <- ascn %>%
        dplyr::filter(cell_id == cell_segments$cell_id[i], chr == cell_segments$chr[i]) %>%
        dplyr::mutate(switch = FALSE) %>%
        dplyr::select(chr, start, end, phase, switch)
      print(dim(ascn_temp))
      bin_switch <- bind_rows(bin_switch, ascn_temp)
    } else{
      ascn_temp <- ascn %>%
        dplyr::filter(cell_id == cell_segments$cell_id[i], chr == cell_segments$chr[i]) %>%
        dplyr::mutate(switch = ifelse(phase == "B", TRUE, FALSE)) %>%
        dplyr::select(chr, start, end, phase, switch)
      bin_switch <- bind_rows(bin_switch, ascn_temp)
      print(dim(ascn_temp))
    }
    print("")
  }
  return(bin_switch)
}

rephase_alleles <- function(ascn){
  bin_switch <- find_switches(ascn, binfilter = 10)
  ascn_switch <- as.data.table(ascn)

  ascn_switch <- ascn_switch[bin_switch, on = .(chr, start, end)] %>%
    .[switch==TRUE, c("Maj", "Min") := .(Min, Maj)] %>%
    .[switch==TRUE, BAF := 1 - BAF] %>%
    .[switch==TRUE, c("alleleA", "alleleB") := .(alleleB, alleleA)] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
                                           1 * ((Min < Maj)) +
                                           2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                                 1 * ((Min < Maj) & (Min != 0)) +
                                                                                 2 * ((Min > Maj) & (Maj != 0)) +
                                                                                 3 * ((Min < Maj) & (Min == 0)) +
                                                                                 4 * ((Min > Maj) & (Maj == 0))]
      ] %>%
    #.[, c("Maj", "Min") := NULL] %>%
    .[order(cell_id, chr, start)]

  ascn_switch <- as.data.frame(ascn_switch)

  ascn_switch2 <- callAlleleSpecificCNHMM(ascn_switch, ncores = 2)
  cl <- umap_clustering(ascn, minPts = max(round(0.05 * length(unique(ascn$cell_id))), 2))
  #ascn <- left_join()

  pdf("~/Downloads/heatmap_2.pdf", width = 25)
  plotHeatmap(ascn, clusters = cl$clustering, plotcol = "state", plottree = F, reorderclusters = T)
  plotHeatmap(ascn, clusters = cl$clustering, plotcol = "state_phase", plottree = F, reorderclusters = T)
  plotHeatmap(ascn_switch2, clusters = cl$clustering, plotcol = "state_phase", plottree = F, reorderclusters = T)
  dev.off()
}



