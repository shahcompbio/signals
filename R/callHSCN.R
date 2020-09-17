#' @export
HaplotypeHMM <- function(n,
                      x,
                      binstates,
                      minor_cn,
                      loherror = 0.02,
                      selftransitionprob = 0.999,
                      eps = 1e-12,
                      likelihood = "binomial",
                      rho = 0.0,
                      Abias = 0.0){

  minor_cn_mat <- t(replicate(length(binstates), minor_cn))
  total_cn_mat <- replicate(length(minor_cn), binstates)

  p <- t(vapply(binstates, function(x) minor_cn / x, FUN.VALUE = numeric(length(minor_cn))))
  p[minor_cn_mat < total_cn_mat / 2] <- p[minor_cn_mat < total_cn_mat / 2] + loherror
  p[minor_cn_mat > total_cn_mat / 2] <- p[minor_cn_mat > total_cn_mat / 2] - loherror

  if (likelihood == "binomial"){
    l <- suppressWarnings(dbinom(x, n, p, log = TRUE))
  } else{
    l <- suppressWarnings(VGAM::dbetabinom(x, n, p, rho = rho, log = TRUE))
    dim(l) <- c(length(x), length(minor_cn))
  }
  if (eps > 0.0){
    l <- matrix(vapply(l, function(x) logspace_addcpp(x, log(eps)), FUN.VALUE = numeric(1)), dim(l)[1], dim(l)[2])
  }
  if (Abias > 0.0){
    lA <- l[minor_cn_mat < total_cn_mat / 2]
    lA <- vapply(lA, function(x) logspace_addcpp(x, log(Abias)), FUN.VALUE = numeric(1))
    l[minor_cn_mat < total_cn_mat / 2] <- lA
  }
  l[is.na(l)] <- log(0.0)
  l[minor_cn_mat > total_cn_mat] <- log(0.0)

  if (selftransitionprob == 0.0){
    tProbs <- matrix(1 / length(minor_cn), length(minor_cn), length(minor_cn))
  } else{
    tProbs <- matrix((1 - selftransitionprob) / (length(minor_cn) - 1), length(minor_cn), length(minor_cn))
    diag(tProbs) <- selftransitionprob
    colnames(tProbs) <- paste0(minor_cn)
    row.names(tProbs) <- paste0(minor_cn)
  }

  states <- paste0(minor_cn)

  res <- viterbi(l, log(tProbs), observations = 1:length(binstates))

  return(list(minorcn = res, l = l))
}

#' @export
assignHaplotypeHMM <- function(CNBAF,
                            minor_cn,
                            eps = 1e-12,
                            loherror = 0.02,
                            selftransitionprob = 0.999,
                            pb = NULL,
                            likelihood = "binomial",
                            rho = 0.0,
                            Abias = 0.0){

  if (!is.null(pb)){
    pb$tick()$print()
  }

  minorcn_res <- c()
  for (mychr in unique(CNBAF$chr)){
    hmmresults <- HaplotypeHMM(n = dplyr::filter(CNBAF, chr == mychr)$totalcounts,
                               x = dplyr::filter(CNBAF, chr == mychr)$alleleB,
                               dplyr::filter(CNBAF, chr == mychr)$state,
                               minor_cn,
                               loherror = loherror,
                               eps = eps,
                               selftransitionprob = selftransitionprob,
                               rho = rho,
                               likelihood = likelihood,
                               Abias = Abias)
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
                              selftransitionprob = 0.999){

  minorcn_res <- c()
  for (mychr in unique(CNBAF$chr)){
    hmmresults <- HaplotypeHMM(n = dplyr::filter(CNBAF, chr == mychr)$totalcounts,
                            x = dplyr::filter(CNBAF, chr == mychr)$alleleB,
                            dplyr::filter(CNBAF, chr == mychr)$state,
                            minor_cn,
                            loherror = loherror,
                            eps = eps,
                            selftransitionprob = selftransitionprob,
                            rho = rho,
                            likelihood = likelihood)
    minorcn_res <- c(minorcn_res, hmmresults$minorcn)
  }

  CNBAF$state_min <- as.numeric(minorcn_res)

  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>% #catch edge cases of 0|1 and 1|0 states
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)]

  CNBAF <- add_states(CNBAF)

  return(list(alleleCN = CNBAF, posterior_prob = hmmresults$posterior_prob, l = hmmresults$l))
}

.callHaplotypeSpecificCN_ <- function(CNBAF,
                                    eps = 1e-12,
                                    loherror = 0.02,
                                    maxCN = 10,
                                    selftransitionprob = 0.999,
                                    progressbar = TRUE,
                                    ncores = 1,
                                    likelihood = "binomial",
                                    rho = 0.0,
                                    Abias = 0.0){

  minor_cn <- seq(0, maxCN, 1)

  if (progressbar == TRUE){
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else{
    pb <- NULL
  }

  #Make sure dataframe is in chromosome position order
  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[order(cell_id, chr, start)]

  if (ncores > 1){
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
                                                    function(cell) assignHaplotypeHMM(dplyr::filter(CNBAF, cell_id == cell),
                                                                                      minor_cn,
                                                                                      eps = eps, loherror = loherror,
                                                                                      selftransitionprob = selftransitionprob,
                                                                                      likelihood = likelihood,
                                                                                      rho = rho,
                                                                                      Abias = Abias,
                                                                                      pb = pb), mc.cores = ncores)) %>%
      .[order(cell_id, chr, start)]
  } else{
    alleleCN <- data.table::rbindlist(lapply(unique(CNBAF$cell_id),
                                        function(cell) assignHaplotypeHMM(dplyr::filter(CNBAF, cell_id == cell), minor_cn,
                                                                       eps = eps, loherror = loherror,
                                                                       likelihood = likelihood,
                                                                       rho = rho,
                                                                       Abias = Abias,
                                                                       selftransitionprob = selftransitionprob,
                                                                       pb = pb))) %>%
      .[order(cell_id, chr, start)]
  }

  alleleCN <- alleleCN %>%
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>% #catch edge cases of 0|1 and 1|0 states
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)]

  alleleCN <- add_states(alleleCN)

  return(as.data.frame(alleleCN))
}

min_cells <- function(haplotypes, minfrachaplotypes = 0.95, mincells = NULL){
  nhaps_vec <- c()
  prop_vec <- c()
  prop <- seq(0.025, 1.0, 0.025)
  mycells <- unique(haplotypes$cell_id)
  ncells <- length(mycells)
  haplotype_counts <- as.data.table(haplotypes) %>%
    .[, list(n = .N, nfrac = .N / ncells), by = c("chr", "start", "end", "hap_label")]
  nhaps <- dim(dplyr::distinct(haplotype_counts, chr, start, hap_label))[1]
  for (i in prop){
    samp_cells <- sample(mycells, round(i * length(mycells)))
    haplotype_counts_temp <- as.data.table(dplyr::filter(haplotypes, cell_id %in% samp_cells)) %>%
      .[, list(n = .N, nfrac = .N / ncells), by = c("chr", "start", "end", "hap_label")]
    nhaps_temp <- dim(dplyr::distinct(haplotype_counts_temp, chr, start, hap_label))[1]
    nhaps_vec <- c(nhaps_vec, nhaps_temp)
    prop_vec <- c(prop_vec, i)
  }

  df <- data.frame(prop = prop_vec, nhaps = nhaps_vec / nhaps) %>%
    dplyr::mutate(ncells = round(prop * length(mycells)))

  ncells_forclustering <- df %>%
    dplyr::filter(nhaps > minfrachaplotypes) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(ncells)

  if (is.null(mincells)){
    mincells <- round(0.01 * length(unique(mycells)))
  }

  ncells_forclustering <- max(ncells_forclustering, mincells)
  return(list(ncells_forclustering = ncells_forclustering, prop = df))
}

#' @export
proportion_imbalance_manual <- function(ascn, haplotypes, cl, field = "copy", phasebyarm = FALSE){
  ncells <- length(unique(ascn$cell_id))

  alleles <- data.table()
  ascn <- as.data.table(dplyr::left_join(ascn, cl$clustering))

  #removed cells that are in the unnassigned group
  ascn <- ascn[clone_id != "0"]

  #calculate proportion of bins in each chromosome or chrarm that exhibit allelic imbalance
  if (phasebyarm) {
    message("Phasing by chromosome arm...")
    ascn$chrarm <- paste0(ascn$chr, coord_to_arm(ascn$chr, ascn$start))
    prop <- ascn[, list(propA = sum(balance) / .N), by = .(chrarm, clone_id)]
    prop <- prop[prop[, .I[which.max(propA)], by=chrarm]$V1]
  } else {
    message("Phasing by chromosome (not arm level)..")
    prop <- ascn[, list(propA = sum(balance) / .N), by = .(chr, clone_id)]
    prop <- prop[prop[, .I[which.max(propA)], by=chr]$V1]
  }

  return(list(prop = prop, cl = cl))
}

#' @export
proportion_imbalance <- function(ascn, haplotypes,
                                 field = "copy",
                                 phasebyarm = FALSE,
                                 minfrachaplotypes = 0.95,
                                 clustering_method = "copy",
                                 overwritemincells = NULL){
  ncells <- length(unique(ascn$cell_id))
  if (is.null(overwritemincells)){
    ncells_for_clustering <- min_cells(haplotypes, minfrachaplotypes = minfrachaplotypes)
    ncells_for_clustering <- ncells_for_clustering$ncells_forclustering
  } else{
    ncells_for_clustering <- overwritemincells
  }
  message(paste0("Using ", ncells_for_clustering, " cells for clustering..."))

  #cluster cells using umap and the "copy" corrected read count value
  if (clustering_method == "copy"){
    cl <- umap_clustering(ascn,
                          minPts = ncells_for_clustering,
                          field = field)
  } else {
    cl <- umap_clustering_breakpoints(ascn,
                          minPts = ncells_for_clustering)
  }
  alleles <- data.table()
  ascn <- as.data.table(dplyr::left_join(ascn, cl$clustering))

  #removed cells that are in the unnassigned group
  ascn <- ascn[clone_id != "0"]

  #calculate proportion of bins in each chromosome or chrarm that exhibit allelic imbalance
  if (phasebyarm) {
    message("Phasing by chromosome arm...")
    ascn$chrarm <- paste0(ascn$chr, coord_to_arm(ascn$chr, ascn$start))
    prop <- ascn[, list(propA = round(sum(balance) / .N, 2), n = sum(balance)), by = .(chrarm, clone_id)]
    prop <- prop[order(n, decreasing = TRUE)] # if there are ties make sure the clone with the largest number of cells gets chosen
    prop <- prop[prop[, .I[which.max(propA)], by=chrarm]$V1]
  } else {
    message("Phasing by chromosome (not arm level)..")
    prop <- ascn[, list(propA = round(sum(balance) / .N, 2), n = sum(balance)), by = .(chr, clone_id)]
    prop <- prop[order(n, decreasing = TRUE)]
    prop <- prop[prop[, .I[which.max(propA)], by=chr]$V1]
  }

  return(list(prop = prop, cl = cl))
}

#' @export
phase_haplotypes_bychr <- function(haplotypes, prop, phasebyarm = FALSE){
  haplotypes <- as.data.table(haplotypes)[as.data.table(prop$cl$clustering), on = "cell_id"]

  if (phasebyarm) {
    haplotypes$chrarm <- paste0(haplotypes$chr, coord_to_arm(haplotypes$chr, haplotypes$start))
    phased_haplotypes <- data.table()
    for (i in 1:nrow(prop$prop)){
      phased_haplotypes_temp <- haplotypes[clone_id == prop$prop[i,"clone_id"] & chrarm == prop$prop[i,"chrarm"]] %>%
        .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
        .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
        .[, c("allele1", "allele0") := NULL]
      phased_haplotypes <- rbind(phased_haplotypes, phased_haplotypes_temp)
    }
  }  else {
    phased_haplotypes <- data.table()
    for (i in 1:nrow(prop$prop)){
      phased_haplotypes_temp <- haplotypes[clone_id == prop$prop[i,"clone_id"]$clone_id & chr == prop$prop[i,"chr"]$chr] %>%
        .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
        .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
        .[, c("allele1", "allele0") := NULL]
      phased_haplotypes <- rbind(phased_haplotypes, phased_haplotypes_temp)
    }
  }

  return(phased_haplotypes)
}

tarones_Z <- function(alleleA, totalcounts){
  m_null <- alleleA
  n <- totalcounts
  p_hat = sum(m_null) / sum(n)
  S = sum( (m_null - n * p_hat)^2 / (p_hat * (1 - p_hat)) )
  Z_score = (S - sum(n)) / sqrt(2 * sum(n * (n - 1)))
  return(Z_score)
}

fitBB <- function(ascn){
  modal_state <- which.max(table(dplyr::filter(ascn, Min > 0)$state_AS_phased))
  bdata <- dplyr::filter(ascn, state_AS_phased == names(modal_state))
  nsample <- min(length(bdata$state), 10^5)
  if(nsample == 10^5){
   bdata <- dplyr::sample_n(bdata, 10^5)
  }
  expBAF <- bdata %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::mutate(e = Min / (Min + Maj)) %>%
    dplyr::pull(e)
  message(paste0('Fitting beta-binomial model to state: ', names(modal_state), "..."))
  fit <- VGAM::vglm(cbind(alleleB, totalcounts-alleleB) ~ 1,
                    VGAM::betabinomial,
                    data = bdata,
                    trace = TRUE)
  message(paste0("Inferred mean: ", round(VGAM::Coef(fit)[["mu"]], 3), ", Expected mean: ", round(expBAF, 3), ", Inferred overdispersion (rho): ", round(VGAM::Coef(fit)[["rho"]], 4)))
  rho <- VGAM::Coef(fit)[["rho"]]
  tz <- tarones_Z(bdata$alleleB, bdata$totalcounts)
  bbfit <- list(fit = fit,
                rho = rho,
                likelihood = "betabinomial",
                expBAF = expBAF,
                state = modal_state,
                taronesZ = tz)
  return(bbfit)
}


#' Call haplotype specific copy number in single cell datasets
#'
#' @param CNbins single cell copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`
#' @param haplotypes single cell haplotypes dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`, `totalcounts`
#' @param eps default 1e-12
#' @param loherror LOH error rate for initial assignment, this is inferred directly from the data in the second pass, default = 0.02
#' @param maxCN maximum copy number to infer allele specific states, default=NULL which will use the maximum state from CNbins
#' @param selftransitionprob probability to stay in the same state in the HMM, default = 0.999, set to 0.0 for an IID model
#' @param progressbar Boolean to display progressbar or not, default = TRUE, will only show if ncores == 1
#' @param ncores Number of cores to use, default = 1
#' @param minfrac Minimum proportion of haplotypes to retain when clustering + phasing, default = 0.8
#' @param likelihood Likelihood model for HMM, default is `binomial`, other option is `betabinomial` or use `auto` and the algorithm will choose the likelihood that best fits the data.
#' @param minbins Minimum number of bins containing both haplotype counts and copy number data for a cell to be included
#' @param minbinschr Minimum number of bins containing both haplotype counts and copy number data per chromosome for a cell to be included
#' @param phased_haplotypes Use this if you want to manually define the haplotypes phasing if for example the default heuristics used by schnapps does not return a good fit.
#' @param clustering_method Method to use to cluster cells for haplotype phasing, default is `copy`, other option is `breakpoints`
#' @param maxloherror Maximum value for LOH error rate
#' @param overwritemincells default NULL
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
#' @examples
#' sim_data <- simulate_data_cohort(clone_num = c(20, 20),
#'        clonal_events = list(list("1" = c(2,0), "5" = c(3,1)),
#'                        list("2" = c(6,3), "3" = c(1,0))),
#'        loherror = 0.02,
#'        coverage = 100)
#'
#' results <- callHaplotypeSpecificCN(sim_data$CNbins, sim_data$haplotypes)
#'
#' @md
#' @export
callHaplotypeSpecificCN <- function(CNbins,
                                    haplotypes,
                                    eps = 1e-12,
                                    loherror = 0.02,
                                    maxCN = NULL,
                                    selftransitionprob = 0.999,
                                    progressbar = TRUE,
                                    ncores = 1,
                                    phasebyarm = FALSE,
                                    minfrac = 0.8,
                                    likelihood = "binomial",
                                    minbins = 100,
                                    minbinschr = 10,
                                    phased_haplotypes = NULL,
                                    clustering_method = "copy",
                                    maxloherror = 0.035,
                                    overwritemincells = NULL) {

  if (!clustering_method %in% c("copy", "breakpoints")){
    stop("Clustering method must be one of copy or breakpoints")
  }

  if (!likelihood %in% c("binomial", "betabinomial", "auto")){
    stop("Likelihood model for HMM emission model must be one of binomial, betabinomial or auto",
         call. = FALSE)
  }

  if (likelihood == "betabinomial" | likelihood == "auto"){
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
           call. = FALSE)
    }
  }

  if (is.null(maxCN)){
    maxCN <- max(CNbins$state)
  }


  cnbaf <- combineBAFCN(haplotypes = haplotypes, CNbins = CNbins)
  ascn <- schnapps:::.callHaplotypeSpecificCN_(cnbaf,
                                    eps = eps,
                                    loherror = loherror,
                                    maxCN = maxCN,
                                    selftransitionprob = selftransitionprob,
                                    progressbar = progressbar,
                                    ncores = ncores)

  infloherror <- ascn %>%
    dplyr::filter(state_phase == "A-LOH") %>%
    dplyr::summarise(err = weighted.mean(x = BAF, w = totalcounts, na.rm = TRUE)) %>%
    dplyr::pull(err)
  infloherror <- min(infloherror, maxloherror) #ensure loh error rate is < maxloherror

  if (likelihood == 'betabinomial' | likelihood == "auto"){
    bbfit <- fitBB(ascn)
    if (bbfit$taronesZ > 5){
      likelihood <- "betabinomial"
      message(paste0("Tarones Z-score: ", round(bbfit$taronesZ, 3), ", using ", likelihood, " model for inference."))
    } else {
      likelihood <- "binomial"
      message(paste0("Tarones Z-score: ", round(bbfit$taronesZ, 3), ", using ", likelihood, " model for inference."))
    }
  } else{
    bbfit <- list(fit = NULL,
                  rho = 0.0,
                  likelihood = "binomial",
                  expBAF = NULL,
                  state = NULL,
                  taronesZ = NULL)
  }

  ascn$balance <- ifelse(ascn$phase == "Balanced", 0, 1)

  if (is.null(phased_haplotypes)){
    p <- proportion_imbalance(ascn,
                              haplotypes,
                              phasebyarm = phasebyarm,
                              minfrac = minfrac,
                              clustering_method = clustering_method,
                              overwritemincells = overwritemincells)
    phased_haplotypes <- phase_haplotypes_bychr(haplotypes = haplotypes,
                                                prop = p,
                                                phasebyarm = phasebyarm)
  } else{
    p <- NULL
  }

  cnbaf <- combineBAFCN(haplotypes = haplotypes,
                        CNbins = CNbins,
                        phased_haplotypes = phased_haplotypes,
                        minbinschr = minbinschr,
                        minbins = minbins)

  hscn_data <- schnapps:::.callHaplotypeSpecificCN_(cnbaf,
                                               eps = eps,
                                               loherror = infloherror,
                                               maxCN = maxCN,
                                               selftransitionprob = selftransitionprob,
                                               progressbar = progressbar,
                                               ncores = ncores,
                                               likelihood = likelihood,
                                               rho = bbfit$rho)

  # Output
  out = list()
  class(out) <- "hscn"

  out[["data"]] <- hscn_data %>% as.data.frame()
  out[["phasing"]] <- p
  out[["loherror"]] <- infloherror
  out[["eps"]] <- eps
  out[["likelihood"]] <- bbfit
  out[["qc_summary"]] <- qc_summary(hscn_data)
  out[["haplotype_phasing"]] <- phased_haplotypes

  out <- fix_assignments(out)

  return(out)
}

#' @export
fix_assignments <- function(hscn){

  if (hscn$likelihood$likelihood == "binomial"){
  hscn_data <- hscn$data %>%
    as.data.table() %>%
    .[, rlid := data.table::rleid(state_AS_phased)] %>%
    .[, alleleAtot := sum(alleleA), by = "rlid"] %>%
    .[, alleleBtot := sum(alleleB), by = "rlid"] %>%
    .[, pMin := Min/state] %>%
    .[, pMin := fifelse(pMin == 0.0, pMin + hscn$loherror, pMin)] %>%
    .[, pMin := fifelse(pMin == 1.0, pMin - hscn$loherror, pMin)] %>%
    .[, pMaj := Maj/state] %>%
    .[, pMaj := fifelse(pMaj == 0.0, pMaj + hscn$loherror, pMaj)] %>%
    .[, pMaj := fifelse(pMaj == 1.0, pMaj - hscn$loherror, pMaj)] %>%
    .[, LLassigned := dbinom(alleleAtot, alleleAtot + alleleBtot, p = pMin)] %>%
    .[, LLother := dbinom(alleleAtot, alleleAtot + alleleBtot, p = pMaj)] %>%
    .[, state_min := fifelse(LLother < LLassigned, Maj, Min)] %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)] %>%
    add_states() %>%
    dplyr::select(-LLassigned, LLother, -state_min, -alleleBtot, -alleleAtot, -pMin, -pMaj)

  } else{
    hscn_data <- hscn$data %>%
      as.data.table() %>%
      .[, rlid := data.table::rleid(state_AS_phased)] %>%
      .[, alleleAtot := sum(alleleA), by = "rlid"] %>%
      .[, alleleBtot := sum(alleleB), by = "rlid"] %>%
      .[, pMin := Min/state] %>%
      .[, pMin := fifelse(pMin == 0.0, pMin + hscn$loherror, pMin)] %>%
      .[, pMin := fifelse(pMin == 1.0, pMin - hscn$loherror, pMin)] %>%
      .[, pMaj := Maj/state] %>%
      .[, pMaj := fifelse(pMaj == 0.0, pMaj + hscn$loherror, pMaj)] %>%
      .[, pMaj := fifelse(pMaj == 1.0, pMaj - hscn$loherror, pMaj)] %>%
      .[, LLassigned := VGAM::dbetabinom(alleleAtot, alleleAtot + alleleBtot, rho = hscn$likelihood$rho, p = pMin)] %>%
      .[, LLother := VGAM::dbetabinom(alleleAtot, alleleAtot + alleleBtot,rho = hscn$likelihood$rho, p = pMaj)] %>%
      .[, state_min := fifelse(LLother < LLassigned, Maj, Min)] %>%
      .[, Maj := state - state_min] %>%
      .[, Min := state_min] %>%
      .[, Min := fifelse(Min < 0, 0, Min)] %>%
      .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
      .[, Min := fifelse(Min > state, state, Min)] %>%
      .[, Maj := fifelse(Maj > state, state, Maj)] %>%
      add_states() %>%
      dplyr::select(-LLassigned, LLother, -state_min, -alleleBtot, -alleleAtot, -pMin, -pMaj)
  }

  hscn[["data"]] <- hscn_data %>% as.data.frame()
  hscn[["qc_summary"]] <- qc_summary(hscn_data)
  return(hscn)
}
