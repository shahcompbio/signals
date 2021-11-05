logspace_add <- function(logx, logy) {
  pmax(logx, logy) + log1p(exp(-abs(logx - logy)))
}

#' @export
alleleHMM <- function(n,
                      x,
                      binstates,
                      minor_cn,
                      loherror = 0.02,
                      selftransitionprob = 0.999,
                      eps = 1e-12,
                      rho = 0.0,
                      likelihood = "binomial") {
  minor_cn_mat <- t(replicate(length(binstates), minor_cn))
  total_cn_mat <- replicate(length(minor_cn), binstates)

  p <- t(vapply(binstates, function(x) minor_cn / x,
    FUN.VALUE = numeric(length(minor_cn))
  ))
  p[minor_cn_mat < total_cn_mat / 2] <-
    p[minor_cn_mat < total_cn_mat / 2] + loherror
  p[minor_cn_mat > total_cn_mat / 2] <-
    p[minor_cn_mat > total_cn_mat / 2] - loherror

  if (likelihood == "binomial") {
    l1log <- suppressWarnings(dbinom(x, n, p, log = T))
    l2log <- suppressWarnings(dbinom(n - x, n, p, log = T))
    l1l2log <- mapply(function(x, y) logspace_addcpp(x, y), l1log, l2log)
    l <- structure(l1l2log, dim = dim(l1log))
  } else {
    l1log <- suppressWarnings(VGAM::dbetabinom(x, n, p, rho = rho, log = T))
    dim(l1log) <- c(length(x), length(minor_cn))
    l2log <- suppressWarnings(VGAM::dbetabinom(n - x, n, p, rho = rho, log = T))
    dim(l2log) <- c(length(x), length(minor_cn))
    l1l2log <- mapply(function(x, y) logspace_addcpp(x, y), l1log, l2log)
    l <- structure(l1l2log, dim = dim(l1log))
  }
  if (eps > 0.0) {
    ltemp <- vapply(l, function(x) logspace_addcpp(x, log(eps)),
      FUN.VALUE = double(1)
    )
    l <- matrix(ltemp, dim(l)[1], dim(l)[2])
  }
  l[is.na(l)] <- log(0.0)
  l[minor_cn_mat > total_cn_mat / 2] <- log(0.0)

  if (selftransitionprob == 0.0) {
    transition_prob <- matrix(
      1 / length(minor_cn),
      length(minor_cn), length(minor_cn)
    )
  } else {
    transition_prob <- matrix(
      (1 - selftransitionprob) / (length(minor_cn) - 1),
      length(minor_cn), length(minor_cn)
    )
    diag(transition_prob) <- selftransitionprob
    colnames(transition_prob) <- paste0(minor_cn)
    row.names(transition_prob) <- paste0(minor_cn)
  }

  res <- viterbi(l, log(transition_prob),
    observations = seq_len(length(binstates))
  )

  return(list(minorcn = res, l = l))
}

#' @export
assignalleleHMM <- function(CNBAF,
                            minor_cn,
                            eps = 1e-12,
                            loherror = 0.02,
                            selftransitionprob = 0.999,
                            pb = NULL,
                            rho = 0.0,
                            likelihood = "binomial") {
  if (!is.null(pb)) {
    pb$tick()$print()
  }

  minorcn_res <- c()
  for (mychr in unique(CNBAF$chr)) {
    hmmresults <- alleleHMM(
      n = dplyr::filter(CNBAF, chr == mychr)$totalcounts,
      x = dplyr::filter(CNBAF, chr == mychr)$alleleB,
      dplyr::filter(CNBAF, chr == mychr)$state,
      minor_cn,
      loherror = loherror,
      eps = eps,
      selftransitionprob = selftransitionprob,
      rho = rho,
      likelihood = likelihood
    )
    minorcn_res <- c(minorcn_res, hmmresults$minorcn)
  }

  CNBAF$state_min <- as.numeric(minorcn_res)

  CNBAF <- data.table::as.data.table(CNBAF)

  return(as.data.frame(CNBAF))
}

#' @export
callalleleHMMcell <- function(CNBAF,
                              minor_cn,
                              eps = 1e-12,
                              loherror = 0.02,
                              selftransitionprob = 0.999) {
  hmmresults <- alleleHMM(
    n = CNBAF$totalcounts,
    x = CNBAF$alleleB,
    CNBAF$state,
    minor_cn,
    loherror = loherror,
    eps = eps,
    selftransitionprob = selftransitionprob
  )

  CNBAF$state_min <- as.numeric(hmmresults$minorcn)

  CNBAF <- data.table::as.data.table(CNBAF)

  return(list(
    alleleCN = CNBAF,
    posterior_prob = hmmresults$posterior_prob,
    l = hmmresults$l
  ))
}

switch_alleles <- function(cn) {
  phase_cn <- cn %>%
    as.data.table() %>%
    .[, switch := data.table::fifelse(
      Min > Maj,
      "switch",
      "stick"
    )] %>%
    .[, alleleA := data.table::fifelse(
      switch == "switch",
      alleleB,
      alleleA
    )] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, switch := data.table::fifelse(
      phase != "Balanced" & BAF > 0.5,
      "switch",
      "stick"
    )] %>%
    .[, alleleA := data.table::fifelse(
      switch == "switch",
      alleleB,
      alleleA
    )] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, c(
      "chr", "start", "end", "cell_id", "state",
      "copy", "alleleA", "alleleB", "totalcounts"
    )] %>%
    .[, BAF := alleleB / totalcounts]
  return(phase_cn)
}

#' Call allele specific copy number in single cell datasets
#'
#' @param CNbins single cell copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`
#' @param haplotypes single cell haplotypes dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`, `totalcounts`
#' @param eps default 1e-12
#' @param loherror LOH error rate for initial assignment, this is inferred directly from the data in the second pass, default = 0.02
#' @param maxCN maximum copy number to infer allele specific states, default=NULL which will use the maximum state from CNbins
#' @param selftransitionprob probability to stay in the same state in the HMM, default = 0.999, set to 0.0 for an IID model
#' @param progressbar Boolean to display progressbar or not, default = TRUE, will only show if ncores == 1
#' @param ncores Number of cores to use, default = 1
#' @param likelihood Likelihood model for HMM, default is `binomial`, other option is `betabinomial` or use `auto` and the algorithm will choose the likelihood that best fits the data.
#' @param minbins Minimum number of bins containing both haplotype counts and copy number data for a cell to be included
#' @param minbinschr Minimum number of bins containing both haplotype counts and copy number data per chromosome for a cell to be included
#' @param maxloherror Maximum value for LOH error rate
#' @param filterhaplotypes filter out haplotypes present in less than X fraction, default is 0.1
#' @param fillmissing For bins with missing counts fill in values based on neighbouring bins
#'
#' @return allele specific copy number object which includes dataframe similar to input with additional columns which include
#'
#' * `Maj` (Major allele copy number)
#' * `Min` (Minor allele copy number)
#' * `state_AS_phased` (phased state of the form Maj|Min )
#' * `state_AS` (mirrored state of the form Maj|Min)
#' * `LOH` (is bin LOH or not)
#' * `state_phase` (state describing which is the dominant allele and whether it is LOH or not)
#' * `state_BAF` (binned discretized BAF value calculated as Min / (Maj + Min))
#' 
#'
#' @details
#' In the allele specific copy number inference Maj is always > Min and state_AS_phased == state_AS
#'
#' @export
callAlleleSpecificCN <- function(CNbins,
                                 haplotypes,
                                 eps = 1e-12,
                                 loherror = 0.02,
                                 maxCN = NULL,
                                 selftransitionprob = 0.95,
                                 progressbar = TRUE,
                                 ncores = 1,
                                 likelihood = "binomial",
                                 minbins = 100,
                                 minbinschr = 10,
                                 maxloherror = 0.03,
                                 filterhaplotypes = 0.1, 
                                 fillmissing = TRUE) {
  if (!likelihood %in% c("binomial", "betabinomial", "auto")) {
    stop("Likelihood model for HMM emission model must
         be one of binomial, betabinomial or auto",
      call. = FALSE
    )
  }

  if (likelihood == "betabinomial" | likelihood == "auto") {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the
           beta-binomial model. Please install it.",
        call. = FALSE
      )
    }
  }

  if (is.null(maxCN)) {
    maxCN <- max(CNbins$state)
  }
  
  if (filterhaplotypes){
    haplotypes <- filter_haplotypes(haplotypes, filterhaplotypes)
  }

  CNBAF <- combineBAFCN(haplotypes = haplotypes, CNbins = CNbins)

  # Make sure dataframe is in chromosome position order
  CNBAF <- CNBAF %>%
    data.table::as.data.table() %>%
    .[order(cell_id, chr, start)]

  message("Initial assignment...")
  hscn <- .callHaplotypeSpecificCN_(CNBAF,
    eps = eps,
    loherror = loherror,
    maxCN = maxCN,
    selftransitionprob = selftransitionprob,
    progressbar = progressbar,
    ncores = ncores
  )

  infloherror <- hscn %>%
    dplyr::filter(state_phase == "A-Hom") %>%
    dplyr::summarise(err = weighted.mean(
      x = BAF,
      w = totalcounts,
      na.rm = TRUE
    )) %>%
    # ensure BAF calculations with low counts don't overwhelm signal
    dplyr::pull(err)
  infloherror <- min(infloherror, maxloherror) # ensure loh error rate is < maxloherror

  if (likelihood == "betabinomial" | likelihood == "auto") {
    bbfit <- fitBB(hscn)
    if (bbfit$taronesZ > 5) {
      likelihood <- "betabinomial"
      message(paste0(
        "Tarones Z-score: ", round(bbfit$taronesZ, 3),
        ", using ", likelihood, " model for inference."
      ))
    } else {
      likelihood <- "binomial"
      message(paste0(
        "Tarones Z-score: ", round(bbfit$taronesZ, 3),
        ", using ", likelihood, " model for inference."
      ))
    }
  } else {
    bbfit <- list(
      fit = NULL,
      rho = 0.0,
      likelihood = "binomial",
      expBAF = NULL,
      state = NULL,
      taronesZ = NULL
    )
  }

  CNBAF <- switch_alleles(hscn)
  minor_cn <- seq(0, maxCN, 1)

  message("Rerun using mirrored BAF...")

  if (progressbar == TRUE) {
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else {
    pb <- NULL
  }

  if (ncores > 1) {
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
      function(cell) {
        assignalleleHMM(dplyr::filter(CNBAF, cell_id == cell),
          minor_cn,
          eps = eps,
          loherror = infloherror,
          selftransitionprob = selftransitionprob,
          likelihood = likelihood,
          rho = bbfit$rho,
          pb = pb
        )
      },
      mc.cores = ncores
    )) %>%
      .[order(cell_id, chr, start)]
  } else {
    alleleCN <- data.table::rbindlist(lapply(
      unique(CNBAF$cell_id),
      function(cell) {
        assignalleleHMM(dplyr::filter(CNBAF, cell_id == cell),
          minor_cn,
          eps = eps,
          loherror = infloherror,
          likelihood = likelihood,
          rho = bbfit$rho,
          selftransitionprob = selftransitionprob,
          pb = pb
        )
      }
    )) %>%
      .[order(cell_id, chr, start)]
  }

  alleleCN <- alleleCN %>%
    .[, Maj1 := state - state_min] %>%
    .[, Min1 := state_min] %>%
    # catch edge cases where Min > Maj:
    .[, Min := fifelse(Min1 > Maj1, Maj1, Min1)] %>%
    .[, Maj := fifelse(Min1 > Maj1, Min1, Maj1)] %>%
    .[, Maj1 := NULL] %>%
    .[, Min1 := NULL] %>%
    # catch edge cases of 0|1 and 1|0 states:
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
      1 * ((Min < Maj)) +
      2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-Hom", "B-Hom")[1 +
      1 * ((Min < Maj) & (Min != 0)) +
      2 * ((Min > Maj) & (Maj != 0)) +
      3 * ((Min < Maj) & (Min == 0)) +
      4 * ((Min > Maj) & (Maj == 0))]] %>%
    .[order(cell_id, chr, start)] %>%
    .[, state_BAF := round((Min / state) / 0.1) * 0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]

  # mirror BAF
  alleleCN <- alleleCN[, switch := data.table::fifelse(
    Min > Maj,
    "switch",
    "stick"
  )] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, distA := abs(BAF - (Min / state))] %>%
    .[, distB := abs(BAF - (Maj / state))] %>%
    .[, switch := fifelse((distB < distA) & (BAF != 0.5), "switch", "stick")] %>%
    .[, alleleA := data.table::fifelse(switch == "switch", alleleB, alleleA)] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, BAF := alleleB / (totalcounts)] %>%
    .[, switch := NULL] %>%
    .[, distA := NULL] %>%
    .[, distB := NULL]
  
  
  if (fillmissing){
    alleleCN <- dplyr::left_join(CNbins, alleleCN)
    alleleCN <- tidyr::fill(alleleCN, c("state_min", "Maj", "Min", "state_phase", "state_AS", 
                                        "state_AS_phased", "LOH", "state_BAF", "phase"), .direction = "downup")
  }

  # Output
  out <- list()
  class(out) <- "ascn"

  out[["data"]] <- as.data.frame(alleleCN)
  out[["loherror"]] <- infloherror
  out[["likelihood"]] <- bbfit
  out[["qc_summary"]] <- qc_summary(alleleCN)
  out[["qc_per_cell"]] <- qc_per_cell(alleleCN)

  return(out)
}



#' Call allele specific copy number in single cell datasets
#'
#' @param hscn hscn object from callHaplotypeSpecificCN
#' @param eps default 1e-12
#' @param maxCN maximum copy number to infer allele specific states, default=NULL which will use the maximum state from CNbins
#' @param selftransitionprob probability to stay in the same state in the HMM, default = 0.999, set to 0.0 for an IID model
#' @param progressbar Boolean to display progressbar or not, default = TRUE, will only show if ncores == 1
#' @param ncores Number of cores to use, default = 1
#' @param fillmissing For bins with missing counts fill in values based on neighbouring bins
#'
#' @return allele specific copy number object which includes dataframe similar to input with additional columns which include
#'
#' * `Maj` (Major allele copy number)
#' * `Min` (Minor allele copy number)
#' * `state_AS_phased` (phased state of the form Maj|Min )
#' * `state_AS` (mirrored state of the form Maj|Min)
#' * `LOH` (is bin LOH or not)
#' * `state_phase` (state describing which is the dominant allele and whether it is LOH or not)
#' * `state_BAF` (binned discretized BAF value calculated as Min / (Maj + Min))
#'
#' @details
#' In the allele specific copy number inference Maj is always > Min and state_AS_phased == state_AS
#'
#' @export
callAlleleSpecificCNfromHSCN <- function(hscn,
                                 eps = 1e-12,
                                 maxCN = NULL,
                                 selftransitionprob = 0.95,
                                 progressbar = TRUE,
                                 ncores = 1,
                                 fillmissing = TRUE) {
  
  likelihood <- hscn$likelihood$likelihood
  
  if (likelihood == "betabinomial" | likelihood == "auto") {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the
           beta-binomial model. Please install it.",
           call. = FALSE
      )
    }
  }
  
  if (is.null(maxCN)) {
    maxCN <- max(hscn$data$state)
  }
  
  CNbins <- hscn$data %>% 
    dplyr::select(cell_id, chr, start, end, state, copy)
  
  infloherror <- hscn$loherror
  
  CNBAF <- switch_alleles(hscn$data)
  minor_cn <- seq(0, maxCN, 1)
  
  message("Rerun using mirrored BAF...")
  
  if (progressbar == TRUE) {
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else {
    pb <- NULL
  }
  
  if (ncores > 1) {
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
                                                         function(cell) {
                                                           assignalleleHMM(dplyr::filter(CNBAF, cell_id == cell),
                                                                           minor_cn,
                                                                           eps = eps,
                                                                           loherror = infloherror,
                                                                           selftransitionprob = selftransitionprob,
                                                                           likelihood = hscn$likelihood$likelihood,
                                                                           rho = hscn$likelihood$rho,
                                                                           pb = pb
                                                           )
                                                         },
                                                         mc.cores = ncores
    )) %>%
      .[order(cell_id, chr, start)]
  } else {
    alleleCN <- data.table::rbindlist(lapply(
      unique(CNBAF$cell_id),
      function(cell) {
        assignalleleHMM(dplyr::filter(CNBAF, cell_id == cell),
                        minor_cn,
                        eps = eps,
                        loherror = infloherror,
                        likelihood = hscn$likelihood$likelihood,
                        rho = hscn$likelihood$rho,
                        selftransitionprob = selftransitionprob,
                        pb = pb
        )
      }
    )) %>%
      .[order(cell_id, chr, start)]
  }
  
  alleleCN <- alleleCN %>%
    .[, Maj1 := state - state_min] %>%
    .[, Min1 := state_min] %>%
    # catch edge cases where Min > Maj:
    .[, Min := fifelse(Min1 > Maj1, Maj1, Min1)] %>%
    .[, Maj := fifelse(Min1 > Maj1, Min1, Maj1)] %>%
    .[, Maj1 := NULL] %>%
    .[, Min1 := NULL] %>%
    # catch edge cases of 0|1 and 1|0 states:
    .[, state_min := Min] %>% 
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
                                           1 * ((Min < Maj)) +
                                           2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-Hom", "B-Hom")[1 +
                                                                                 1 * ((Min < Maj) & (Min != 0)) +
                                                                                 2 * ((Min > Maj) & (Maj != 0)) +
                                                                                 3 * ((Min < Maj) & (Min == 0)) +
                                                                                 4 * ((Min > Maj) & (Maj == 0))]] %>%
    .[order(cell_id, chr, start)] %>%
    .[, state_BAF := round((Min / state) / 0.1) * 0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]
  
  # mirror BAF
  alleleCN <- alleleCN[, switch := data.table::fifelse(
    Min > Maj,
    "switch",
    "stick"
  )] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, origBAF := BAF] %>% 
    .[, distA := abs(BAF - (Min / state))] %>%
    .[, distB := abs(BAF - (Maj / state))] %>%
    .[, switch := fifelse((distB < distA) & (BAF != 0.5), "switch", "stick")] %>%
    .[, alleleA := data.table::fifelse(switch == "switch", alleleB, alleleA)] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, BAF := alleleB / (totalcounts)] %>%
    .[, switch := NULL] %>%
    .[, distA := NULL] %>%
    .[, distB := NULL]
  
  if (fillmissing){
    alleleCN <- dplyr::left_join(CNbins, alleleCN)
    alleleCN <- tidyr::fill(alleleCN, c("state_min", "Maj", "Min", "state_phase", "state_AS", 
                                          "state_AS_phased", "LOH", "state_BAF", "phase"), .direction = "downup")
  }
  
  # Output
  out <- list()
  class(out) <- "ascn"
  
  out[["data"]] <- as.data.frame(alleleCN)
  out[["loherror"]] <- infloherror
  out[["likelihood"]] <- hscn$bbfit
  out[["qc_summary"]] <- qc_summary(alleleCN)
  out[["qc_per_cell"]] <- qc_per_cell(alleleCN)
  
  return(out)
}
