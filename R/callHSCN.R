#' @export
callHaplotypeSpecificCN <- function(CNBAF, maxCN = 11){
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
HaplotypeHMM <- function(n,
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
assignHaplotypeHMM <- function(CNBAF,
                            minor_cn,
                            eps = 1e-12,
                            loherror = 0.01,
                            selftransitionprob = 0.999,
                            pb = NULL){

  if (!is.null(pb)){
    pb$tick()$print()
  }

  hmmresults <- HaplotypeHMM(n = CNBAF$totalcounts,
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

  hmmresults <- HaplotypeHMM(n = CNBAF$totalcounts,
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

.callHaplotypeSpecificCN_ <- function(CNBAF,
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
                                                    function(cell) assignHaplotypeHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                                   eps = eps, loherror = loherror,
                                                                                   selftransitionprob = selftransitionprob,
                                                                                   pb = pb), mc.cores = ncores)) %>%
      .[order(cell_id, chr, start)]
  } else{
    alleleCN <- data.table::rbindlist(lapply(unique(CNBAF$cell_id),
                                        function(cell) assignHaplotypeHMM(CNBAF %>% dplyr::filter(cell_id == cell), minor_cn,
                                                                       eps = eps, loherror = loherror,
                                                                       selftransitionprob = selftransitionprob,
                                                                       pb = pb))) %>%
      .[order(cell_id, chr, start)]
  }

  return(as.data.frame(alleleCN))
}

min_cells <- function(haplotypes, minfrachaplotypes = 0.95){
  nhaps_vec <- c()
  prop_vec <- c()
  prop <- seq(0.05, 1.0, 0.01)
  mycells <- unique(haplotypes$cell_id)
  ncells <- length(mycells)
  haplotype_counts <- as.data.table(haplotypes) %>%
    .[, list(n = .N, nfrac = .N / ncells), by = c("chr", "start", "end", "hap_label")]
  nhaps <- dim(dplyr::distinct(haplotype_counts, chr, start, hap_label))[1]
  for (i in prop){
    samp_cells <- sample(mycells, round(i * length(mycells)))
    haplotype_counts_temp <- as.data.table(haplotypes %>% dplyr::filter(cell_id %in% samp_cells)) %>%
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

  ncells_forclustering <- max(ncells_forclustering, round(0.05 * length(unique(mycells))))
  return(ncells_forclustering)
}

#' @export
proportion_imbalance <- function(ascn, haplotypes, field = "copy", phasebyarm = FALSE, minfrac = 0.05){
  ncells <- length(unique(ascn$cell_id))
  ncells_for_clustering <- min_cells(haplotypes)
  message(paste0("Using ", ncells_for_clustering, " cells for clustering..."))

  #cluster cells using umap and the "copy" corrected read count value
  cl <- umap_clustering(ascn,
                        minPts = ncells_for_clustering,
                        field = field)
  alleles <- data.table()
  ascn <- as.data.table(dplyr::left_join(ascn, cl$clustering))

  #removed cells that are in the unnassigned group
  ascn <- ascn[clone_id != "0"]

  #calculate proportion of bins in each chromosome or chrarm that exhibit allelic imbalance
  if (phasebyarm) {
    message("Phasing by chromosome arm...")
    ascn$chrarm <- paste0(ascn$chr, coord_to_arm(ascn$chr, ascn$start))
    prop <- ascn %>%
      .[, list(propA = sum(balance) / .N), by = .(chrarm, clone_id)]
    prop <- prop[prop[, .I[which.max(propA)], by=chrarm]$V1]
  } else {
    message("Phasing by chromosome (not arm level)..")
    prop <- ascn %>%
      .[, list(propA = sum(balance) / .N), by = .(chr, clone_id)]
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


#' @export
callHaplotypeSpecificCN <- function(CNbins,
                                    haplotypes,
                                      eps = 1e-12,
                                      loherror = 0.03,
                                      maxCN = 12,
                                      selftransitionprob = 0.999,
                                      progressbar = TRUE,
                                      ncores = 1,
                                      phasebyarm = FALSE,
                                      minfrac = 0.05) {

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
    dplyr::summarise(err = mean(BAF)) %>%
    dplyr::pull(err)

  ascn$balance <- ifelse(ascn$phase == "Balanced", 0, 1)

  p <- proportion_imbalance(ascn, haplotypes, phasebyarm = phasebyarm, minfrac = minfrac)
  phased_haplotypes <- phase_haplotypes_bychr(haplotypes = haplotypes,
                                              prop = p,
                                              phasebyarm = phasebyarm)
  cnbaf <- combineBAFCN(haplotypes = haplotypes,
                        CNbins = CNbins,
                        phased_haplotypes = phased_haplotypes)

  hscn <- schnapps:::.callHaplotypeSpecificCN_(cnbaf,
                                               eps = eps,
                                               loherror = infloherror,
                                               maxCN = maxCN,
                                               selftransitionprob = selftransitionprob,
                                               progressbar = progressbar,
                                               ncores = ncores)

  # Output
  out = list()
  class(out) <- "hscn"

  out[["data"]] <- hscn %>% dplyr::filter(state_min > -1) # catch weird bug with bin = -1
  out[["phasing"]] <- p
  out[["loherror"]] <- infloherror
  out[["eps"]] <- eps

  return(out)
}
