#' @export
combineBAFCN <- function(haplotypes, CNbins, binsize = 5e6, filtern = 0){

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
  CNbins <- CNbins[haplotypes, on = c("chr", "start", "end", "cell_id")]

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
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                               1 * ((Min < Maj) & (Min != 0)) +
                                                                               2 * ((Min > Maj) & (Maj != 0)) +
                                                                               3 * ((Min < Maj) & (Min == 0)) +
                                                                               4 * ((Min > Maj) & (Maj == 0))]
    ] %>%
  .[, c("Maj", "Min") := NULL] %>%
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

  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
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
