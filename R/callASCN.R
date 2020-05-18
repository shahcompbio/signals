#' @export
alleleHMM <- function(n,
                      x,
                      binstates,
                      minor_cn,
                      loherror = 0.01,
                      selftransitionprob = 0.999,
                      eps = 1e-12,
                      rho = 0.0,
                      likelihood = "binomial"){

  minor_cn_mat <- t(replicate(length(binstates), minor_cn))
  total_cn_mat <- replicate(length(minor_cn), binstates)

  p <- t(vapply(binstates, function(x) minor_cn / x, FUN.VALUE = numeric(length(minor_cn))))
  p[, 1] <- loherror
  p[minor_cn_mat == total_cn_mat] <- 1 - loherror

  if (likelihood == "binomial"){
    l1log <- suppressWarnings(dbinom(x, n, p, log = T))
    l2log <- suppressWarnings(dbinom(n - x, n, p, log = T))
    l1l2log <- mapply(function(x,y) matrixStats::logSumExp(c(x,y)), l1log, l2log)
    l <- structure(l1l2log, dim = dim(l1log))
  } else {
    l1log <- suppressWarnings(VGAM::dbetabinom(x, n, p, rho = rho, log = T))
    dim(l1log) <- c(length(x), length(minor_cn))
    l2log <- suppressWarnings(VGAM::dbetabinom(n - x, n, p, rho = rho, log = T))
    dim(l2log) <- c(length(x), length(minor_cn))
    l1l2log <- mapply(function(x,y) matrixStats::logSumExp(c(x,y)), l1log, l2log)
    l <- structure(l1l2log, dim = dim(l1log))
  }
  if (eps > 0.0){
    l <- matrix(vapply(l, function(x) matrixStats::logSumExp(c(x, log(eps))), FUN.VALUE = numeric(1)), dim(l)[1], dim(l)[2])
  }
  l[is.na(l)] <- log(0.0)
  l[minor_cn_mat > total_cn_mat/2] <- log(0.0)

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
                            pb = NULL,
                            rho = 0.0,
                            likelihood = "binomial"){

  if (!is.null(pb)){
    pb$tick()$print()
  }

  hmmresults <- alleleHMM(n = CNBAF$totalcounts,
                          x = CNBAF$alleleB,
                          CNBAF$state,
                          minor_cn,
                          loherror = loherror,
                          eps = eps,
                          selftransitionprob = selftransitionprob,
                          rho = rho,
                          likelihood = likelihood)

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
  .[order(cell_id, chr, start)] %>%
  .[, state_BAF := round((Min / state)/0.1)*0.1] %>%
  .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]

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
    .[order(cell_id, chr, start)] %>%
    .[, state_BAF := round((Min / state)/0.1)*0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)]

  return(list(alleleCN = CNBAF, posterior_prob = hmmresults$posterior_prob, l = hmmresults$l))
}

switch_alleles <- function(cn, pval = 0.05){
  phase_cn <- cn %>%
    as.data.table() %>%
    .[, switch := data.table::fifelse(Min > Maj, "switch", "stick")] %>%
    .[, alleleA := data.table::fifelse(switch == "switch", alleleB, alleleA)] %>%
    .[, alleleB := totalcounts - alleleA] %>%
    .[, c("chr", "start", "end", "cell_id", "state", "copy","alleleA", "alleleB", "totalcounts")] %>%
    .[, BAF := alleleB / totalcounts]
  return(phase_cn)
}

#' @export
callAlleleSpecificCN <- function(CNbins,
                                 haplotypes,
                                 eps = 1e-12,
                                 loherror = 0.03,
                                 maxCN = 12,
                                 selftransitionprob = 0.999,
                                 progressbar = TRUE,
                                 ncores = 1,
                                 likelihood = "binomial"){

  if (!likelihood %in% c("binomial", "betabinomial")){
    stop("Likelihood model for HMM emission model must be one of binomial and beta-binomial",
         call. = FALSE)
  }

  if (likelihood == "betabinomial"){
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
           call. = FALSE)
    }
  }

  CNBAF <- combineBAFCN(haplotypes = haplotypes, CNbins = CNbins)

  #Make sure dataframe is in chromosome position order
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
                                  ncores = ncores)

  infloherror <- hscn %>%
    dplyr::filter(state_phase == "A-LOH") %>%
    dplyr::summarise(err = mean(BAF, na.rm = T)) %>%
    dplyr::pull(err)

  if (likelihood == 'betabinomial'){
    bbfit <- fitBB(ascn)
  } else{
    bbfit <- list(fit = NULL,
                  rho = 0.0,
                  likelihood = "binomial",
                  expBAF = NULL,
                  state = NULL,
                  taronesZ = NULL)
  }

  CNBAF <- switch_alleles(hscn)
  minor_cn <- seq(0, maxCN, 1)

  message("Rerun using mirrored BAF...")

  if (progressbar == TRUE){
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else{
    pb <- NULL
  }

  if (ncores > 1){
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
                                                    function(cell) assignalleleHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                                   eps = eps,
                                                                                   loherror = infloherror,
                                                                                   selftransitionprob = selftransitionprob,
                                                                                   likelihood = likelihood,
                                                                                   rho = bbfit$rho,
                                                                                   pb = pb), mc.cores = ncores)) %>%
      .[order(cell_id, chr, start)]
  } else{
    alleleCN <- data.table::rbindlist(lapply(unique(CNBAF$cell_id),
                                        function(cell) assignalleleHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                       eps = eps,
                                                                       loherror = infloherror,
                                                                       likelihood = likelihood,
                                                                       rho = bbfit$rho,
                                                                       selftransitionprob = selftransitionprob,
                                                                       pb = pb))) %>%
      .[order(cell_id, chr, start)]
  }

  # Output
  out = list()
  class(out) <- "ascn"

  out[["data"]] <- as.data.frame(alleleCN)
  out[["loherror"]] <- infloherror
  out[["likelihood"]] <- bbfit

  return(out)
}
