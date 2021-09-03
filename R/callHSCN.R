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
                         Abias = 0.0,
                         viterbiver = "cpp") {
  
  #hack to avoid 0/0 numerical errors
  binstates[binstates == 0] <- 1
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
    l <- suppressWarnings(dbinom(x, n, p, log = TRUE))
  } else {
    l <- suppressWarnings(VGAM::dbetabinom(x, n, p, rho = rho, log = TRUE))
    dim(l) <- c(length(x), length(minor_cn))
  }
  if (eps > 0.0) {
    l <- matrix(vapply(l, function(x) logspace_addcpp(x, log(eps)),
                       FUN.VALUE = numeric(1)
    ), dim(l)[1], dim(l)[2])
  }
  if (Abias > 0.0) {
    lA <- l[minor_cn_mat < total_cn_mat / 2]
    lA <- vapply(lA, function(x) logspace_addcpp(x, log(Abias)),
                 FUN.VALUE = numeric(1)
    )
    l[minor_cn_mat < total_cn_mat / 2] <- lA
  }
  l[is.na(l)] <- log(0.0)
  l[minor_cn_mat > total_cn_mat] <- log(0.0)
  
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
  
  if (viterbiver == "R") {
    res <- viterbiR(l, log(transition_prob),
                    observations = seq_len(length(binstates))
    )
  }
  
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
                               Abias = 0.0,
                               viterbiver = "cpp") {
  if (!is.null(pb)) {
    pb$tick()$print()
  }
  
  minorcn_res <- c()
  for (mychr in unique(CNBAF$chr)) {
    hmmresults <- HaplotypeHMM(
      n = dplyr::filter(CNBAF, chr == mychr)$totalcounts,
      x = dplyr::filter(CNBAF, chr == mychr)$alleleB,
      dplyr::filter(CNBAF, chr == mychr)$state,
      minor_cn,
      loherror = loherror,
      eps = eps,
      selftransitionprob = selftransitionprob,
      rho = rho,
      likelihood = likelihood,
      Abias = Abias,
      viterbiver = viterbiver
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
                              selftransitionprob = 0.99) {
  minorcn_res <- c()
  for (mychr in unique(CNBAF$chr)) {
    hmmresults <- HaplotypeHMM(
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
  
  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>%
    # catch edge cases of 0|1 and 1|0 states
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)]
  
  CNBAF <- add_states(CNBAF)
  
  return(list(
    alleleCN = CNBAF,
    posterior_prob = hmmresults$posterior_prob,
    l = hmmresults$l
  ))
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
                                      Abias = 0.0,
                                      viterbiver = "cpp") {
  
  chr <- start <- cell_id <- state_min <- Maj <- Min <- state <- state_AS <- NULL
  state_AS_phased <- LOH <- phase <- state_BAF <- state_phase <- NULL
  
  minor_cn <- seq(0, maxCN, 1)
  
  if (progressbar == TRUE) {
    pb <- dplyr::progress_estimated(length(unique(CNBAF$cell_id)), min_time = 1)
  } else {
    pb <- NULL
  }
  
  # Make sure dataframe is in chromosome position order
  CNBAF <- data.table::as.data.table(CNBAF) %>%
    .[order(cell_id, chr, start)]
  
  if (ncores > 1) {
    alleleCN <- data.table::rbindlist(parallel::mclapply(unique(CNBAF$cell_id),
                                                         function(cell) {
                                                           assignHaplotypeHMM(dplyr::filter(CNBAF, cell_id == cell),
                                                                              minor_cn,
                                                                              eps = eps,
                                                                              loherror = loherror,
                                                                              selftransitionprob = selftransitionprob,
                                                                              likelihood = likelihood,
                                                                              rho = rho,
                                                                              Abias = Abias,
                                                                              pb = pb,
                                                                              viterbiver = viterbiver
                                                           )
                                                         },
                                                         mc.cores = ncores
    )) %>%
      .[order(cell_id, chr, start)]
  } else {
    alleleCN <- data.table::rbindlist(lapply(
      unique(CNBAF$cell_id),
      function(cell) {
        assignHaplotypeHMM(dplyr::filter(CNBAF, cell_id == cell),
                           minor_cn,
                           eps = eps,
                           loherror = loherror,
                           likelihood = likelihood,
                           rho = rho,
                           Abias = Abias,
                           selftransitionprob = selftransitionprob,
                           pb = pb,
                           viterbiver = viterbiver
        )
      }
    )) %>%
      .[order(cell_id, chr, start)]
  }
  
  alleleCN <- alleleCN %>%
    .[, state_min := fifelse(state_min < 0, 0, state_min)] %>%
    # catch edge cases of 0|1 and 1|0 states
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, Min := fifelse(Min < 0, 0, Min)] %>%
    .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
    .[, Min := fifelse(Min > state, state, Min)] %>%
    .[, Maj := fifelse(Maj > state, state, Maj)]
  
  alleleCN <- add_states(alleleCN)
  
  return(as.data.frame(alleleCN))
}


min_cells <- function(haplotypes, minfrachaplotypes = 0.95, mincells = 4, samplen = 5) {
  chr <- hap_label <- NULL
  nhaps_vec <- c()
  prop_vec <- c()
  prop <- c(0.005, 0.01, 0.02, 0.03, 0.04, 0.05, seq(0.1, 1.0, 0.1))
  mycells <- unique(haplotypes$cell_id)

  prop <- c(1 / length(mycells), 2 / length(mycells), 3 / length(mycells), 4 / length(mycells),
            5 / length(mycells), 6 / length(mycells), 7 / length(mycells), 
            8 / length(mycells), 9 / length(mycells), 10 / length(mycells), 
            20 / length(mycells), 50 / length(mycells), prop)
  prop <- unique(sort(prop))
  prop <- prop[prop < 1.0]
  ncells <- length(mycells)
  haplotype_counts <- as.data.table(haplotypes) %>%
    .[, list(n = .N, nfrac = .N / ncells),
      by = c("chr", "start", "end", "hap_label")
      ]
  nhaps <- dim(dplyr::distinct(haplotype_counts, chr, start, hap_label))[1]
  
  haplotype <- as.data.table(haplotypes)
  
  for (i in prop) {
    for (j in 1:samplen){
      if (round(i * length(mycells)) == 0){
        next
      }
      #sample cells
      samp_cells <- sample(mycells, round(i * length(mycells)))
      #filter haplotypes down to sampled cells
      haplotypes_temp <- haplotypes[cell_id %in% samp_cells]
      #count the number of haplotype blocks that have been retained
      nhaps_temp <- dim(dplyr::distinct(haplotypes_temp, chr, start, hap_label))[1]
      nhaps_vec <- c(nhaps_vec, nhaps_temp)
      prop_vec <- c(prop_vec, i)
    }
  }
  
  df <- data.frame(prop = prop_vec, nhaps = nhaps_vec / nhaps) %>%
    dplyr::mutate(ncells = round(prop * length(mycells))) %>% 
    dplyr::group_by(prop, ncells) %>% 

    dplyr::summarise(min_nhaps = min(nhaps), nhaps = mean(nhaps)) %>% 
    dplyr::ungroup(.) %>% 
    as.data.frame(.)

  ncells_forclustering <- df %>%
    dplyr::filter(nhaps >= minfrachaplotypes) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(ncells)
  
  if (is.null(mincells)) {
    mincells <- round(0.01 * length(unique(mycells)))
    mincells <- max(mincells, 5)
  }
  
  ncells_forclustering <- max(ncells_forclustering, mincells)
  return(list(ncells_forclustering = ncells_forclustering, prop = df))
}

#' @export
proportion_imbalance_manual <- function(ascn,
                                        haplotypes,
                                        cl,
                                        field = "copy",
                                        phasebyarm = FALSE) {
  clone_id <- chrarm <- propA <- chr <- NULL
  
  ascn <- as.data.table(dplyr::left_join(ascn, cl$clustering))
  
  # removed cells that are in the unnassigned group
  ascn <- ascn[clone_id != "0"]
  
  # calculate proportion of bins in each chromosome or chrarm that exhibit allelic imbalance
  if (phasebyarm) {
    message("Phasing by chromosome arm...")
    ascn$chrarm <- paste0(ascn$chr, coord_to_arm(ascn$chr, ascn$start))
    prop <- ascn[, list(propA = sum(balance) / .N), by = .(chrarm, clone_id)]
    prop <- prop[prop[, .I[which.max(propA)], by = chrarm]$V1]
  } else {
    message("Phasing by chromosome (not arm level)..")
    prop <- ascn[, list(propA = sum(balance) / .N), by = .(chr, clone_id)]
    prop <- prop[prop[, .I[which.max(propA)], by = chr]$V1]
  }
  
  return(list(prop = prop, cl = cl))
}

get_cells_per_chr_global <- function(ascn,
                                     haplotypes,
                                     ncells_for_clustering,
                                     field = "copy",
                                     clustering_method = "copy",
                                     phasebyarm = FALSE) {
  
  # cluster cells using umap and the "copy" corrected read count value
  if (clustering_method == "copy") {
    cl <- umap_clustering(ascn,
                          minPts = ncells_for_clustering,
                          field = field
    )
  } else {
    cl <- umap_clustering_breakpoints(ascn,
                                      minPts = ncells_for_clustering
    )
  }
  ascn <- as.data.table(dplyr::left_join(ascn, cl$clustering))
  
  # removed cells that are in the unnassigned group
  ascn <- ascn[clone_id != "0"]
  
  # calculate proportion of bins in each chromosome or chrarm that exhibit allelic imbalance
  if (phasebyarm) {
    message("Phasing by chromosome arm...")
    ascn$chrarm <- paste0(ascn$chr, coord_to_arm(ascn$chr, ascn$start))
    prop <- ascn[, list(propA = round(sum(balance) / .N, 2), n = sum(balance)),
                 by = .(chrarm, clone_id)
                 ]
    prop <- prop[order(n, decreasing = TRUE)] # if there are ties make sure the clone with the largest number of cells gets chosen
    prop <- prop[prop[, .I[which.max(propA)], by = chrarm]$V1]
    chrlist <- prop_to_list(as.data.table(haplotypes)[as.data.table(cl$clustering), on = "cell_id"],
                            prop,
                            phasebyarm = phasebyarm
    )
  } else {
    message("Phasing by chromosome (not arm level)..")
    prop <- ascn[, list(propA = round(sum(balance) / .N, 2), n = sum(balance)),
                 by = .(chr, clone_id)
                 ]
    prop <- prop[order(n, decreasing = TRUE)]
    prop <- prop[prop[, .I[which.max(propA)], by = chr]$V1]
    chrlist <- prop_to_list(as.data.table(haplotypes)[as.data.table(cl$clustering), on = "cell_id"],
                            prop,
                            phasebyarm = phasebyarm
    )
  }
  
  return(chrlist)
}

get_cells_per_chr_local <- function(ascn,
                                    haplotypes,
                                    ncells_for_clustering,
                                    field = "state_BAF",
                                    phasebyarm = FALSE) {
  
  # cluster cells per chromosome
  
  chrlist <- list()
  for (mychr in unique(ascn$chr)) {
    message(paste0("Clustering chromosome ", mychr))
    ascn_chr <- as.data.table(ascn)[chr == mychr]
    cl <- umap_clustering(ascn_chr,
                          n_neighbors = 20,
                          min_dist = 0.001,
                          minPts = ncells_for_clustering,
                          field = "state_BAF",
                          umapmetric = "euclidean"
    )
    prop <- ascn_chr[as.data.table(cl$clustering), on = "cell_id"] %>%
      .[, list(
        propA = round(sum(balance) / .N, 2),
        n = sum(balance),
        propModestate = sum(state == Mode(state)) / .N,
        propLOH = sum(LOH == "LOH") / .N,
        ncells = length(unique(cell_id))
      ), by = .(chr, cell_id, clone_id)] %>%
      .[, list(
        propA = median(propA),
        n = median(n),
        propModestate = median(propModestate),
        propLOH = median(propLOH),
        ncells = median(ncells)
      ), by = .(chr, clone_id)]
    prop <- prop[order(propA, propModestate, ncells, propLOH, n, decreasing = TRUE)]
    prop <- prop[prop[, .I[which.max(propA)], by = chr]$V1]
    cells <- dplyr::filter(cl$clustering, clone_id == prop$clone_id[1]) %>%
      dplyr::pull(cell_id)
    chrlist[[mychr]] <- cells
  }
  
  return(chrlist)
}

#' @export
proportion_imbalance <- function(ascn,
                                 haplotypes,
                                 field = "copy",
                                 phasebyarm = FALSE,
                                 clustering_method = "copy",
                                 minfrachaplotypes = 0.95,
                                 mincells = 5,
                                 overwritemincells = NULL,
                                 cluster_per_chr = TRUE) {
  ncells <- length(unique(ascn$cell_id))
  if (is.null(overwritemincells)) {
    ncells_for_clustering <- min_cells(haplotypes,
                                       mincells = mincells,
                                       minfrachaplotypes = minfrachaplotypes
    )
    propdf <- ncells_for_clustering$prop
    ncells_for_clustering <- ncells_for_clustering$ncells_forclustering
  } else {
    ncells_for_clustering <- overwritemincells
    propdf <- NULL
  }
  message(paste0("Using ", ncells_for_clustering, " cells for clustering..."))
  
  if (cluster_per_chr) {
    chrlist <- get_cells_per_chr_local(ascn,
                                       haplotypes,
                                       ncells_for_clustering,
                                       field = field
    )
  } else {
    chrlist <- get_cells_per_chr_global(ascn,
                                        haplotypes,
                                        ncells_for_clustering,
                                        field = field,
                                        clustering_method = clustering_method
    )
  }
  return(list(chrlist = chrlist, propdf = propdf))
}

prop_to_list <- function(haplotypes, prop, phasebyarm = FALSE) {
  if (phasebyarm == TRUE) {
    haplotypes_keep <- dplyr::inner_join(haplotypes, prop,
                                         by = c("chrarm", "clone_id")
    )
    chrlist <- split(haplotypes_keep$cell_id, haplotypes_keep$chrarm)
  } else {
    haplotypes_keep <- dplyr::inner_join(haplotypes, prop,
                                         by = c("chr", "clone_id")
    )
    chrlist <- split(haplotypes_keep$cell_id, haplotypes_keep$chr)
  }
  return(chrlist)
}

#' @export
phase_haplotypes_bychr <- function(haplotypes, chrlist, phasebyarm = FALSE) {
  haplotypes <- as.data.table(haplotypes)
  
  if (phasebyarm) {
    haplotypes$chrarm <- paste0(haplotypes$chr, coord_to_arm(haplotypes$chr, haplotypes$start))
    phased_haplotypes <- data.table()
    for (i in names(chrlist)) {
      phased_haplotypes_temp <- haplotypes[cell_id %in% chrlist[[i]] & chrarm == i] %>%
        .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
        .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
        .[, c("allele1", "allele0") := NULL]
      phased_haplotypes <- rbind(phased_haplotypes, phased_haplotypes_temp)
    }
  } else {
    phased_haplotypes <- data.table()
    for (i in names(chrlist)) {
      phased_haplotypes_temp <- haplotypes[cell_id %in% chrlist[[i]] & chr == i] %>%
        .[, lapply(.SD, sum), by = .(chr, start, end, hap_label), .SDcols = c("allele1", "allele0")] %>%
        .[, phase := ifelse(allele0 < allele1, "allele0", "allele1")] %>%
        .[, c("allele1", "allele0") := NULL]
      phased_haplotypes <- rbind(phased_haplotypes, phased_haplotypes_temp)
    }
  }
  
  return(phased_haplotypes)
}

tarones_Z <- function(alleleA, totalcounts) {
  m_null <- alleleA
  n <- totalcounts
  p_hat <- sum(m_null) / sum(n)
  S <- sum((m_null - n * p_hat)^2 / (p_hat * (1 - p_hat)))
  Z_score <- (S - sum(n)) / sqrt(2 * sum(n * (n - 1)))
  return(Z_score)
}

fitBB <- function(ascn) {
  modal_state <- which.max(table(dplyr::filter(ascn, Min > 0)$state_AS_phased))
  bdata <- dplyr::filter(ascn, state_AS_phased == names(modal_state))
  nsample <- min(length(bdata$state), 10^5)
  if (nsample == 10^5) {
    bdata <- dplyr::sample_n(bdata, 10^5)
  }
  expBAF <- bdata %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::mutate(e = Min / (Min + Maj)) %>%
    dplyr::pull(e)
  message(paste0("Fitting beta-binomial model to state: ", names(modal_state), "..."))
  fit <- VGAM::vglm(cbind(alleleB, totalcounts - alleleB) ~ 1,
                    VGAM::betabinomial,
                    data = bdata,
                    trace = TRUE
  )
  message(paste0("Inferred mean: ", round(VGAM::Coef(fit)[["mu"]], 3), ", Expected mean: ", round(expBAF, 3), ", Inferred overdispersion (rho): ", round(VGAM::Coef(fit)[["rho"]], 4)))
  rho <- VGAM::Coef(fit)[["rho"]]
  tz <- tarones_Z(bdata$alleleB, bdata$totalcounts)
  bbfit <- list(
    fit = fit,
    rho = rho,
    likelihood = "betabinomial",
    expBAF = expBAF,
    state = modal_state,
    taronesZ = tz
  )
  return(bbfit)
}

filter_haplotypes <- function(haplotypes, fraction){
  haplotypes <- as.data.table(haplotypes)
  nhaps <- dim(haplotypes)[1]
  total_cells <- length(unique(haplotypes$cell_id))
  message(paste0("Filtering out haplotypes present < ", fraction * 100, "% of cells..."))
  haplotypes <- haplotypes %>% 
    .[, ncells := length(unique(cell_id)), by = .(chr, start, end, hap_label)] %>% 
    .[, fcells := ncells / total_cells] %>% 
    .[fcells > fraction]
  haplotypes <- haplotypes[, c("fcells", "ncells") := NULL]
  fracretained <- round(dim(haplotypes)[1] / nhaps, 3)
  message(paste0("Fraction of haplotypes retained after filtering = ", fracretained))
  return(haplotypes)
}


#' Call haplotype specific copy number in single cell datasets
#'
#' @param CNbins single cell copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`
#' @param haplotypes single cell haplotypes dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `hap_label`, `allele1`, `allele0`, `totalcounts`
#' @param eps default 1e-12
#' @param loherror LOH error rate for initial assignment, this is inferred directly from the data in the second pass, default = 0.02
#' @param maxCN maximum copy number to infer allele specific states, default=NULL which will use the maximum state from CNbins
#' @param selftransitionprob probability to stay in the same state in the HMM, default = 0.95, set to 0.0 for an IID model
#' @param progressbar Boolean to display progressbar or not, default = TRUE, will only show if ncores == 1
#' @param ncores Number of cores to use, default = 1
#' @param minfrac Minimum proportion of haplotypes to retain when clustering + phasing, default = 0.8
#' @param likelihood Likelihood model for HMM, default is `binomial`, other option is `betabinomial` or use `auto` and the algorithm will choose the likelihood that best fits the data.
#' @param minbins Minimum number of bins containing both haplotype counts and copy number data for a cell to be included
#' @param minbinschr Minimum number of bins containing both haplotype counts and copy number data per chromosome for a cell to be included
#' @param phased_haplotypes Use this if you want to manually define the haplotypes phasing if for example the default heuristics used by schnapps does not return a good fit.
#' @param clustering_method Method to use to cluster cells for haplotype phasing, default is `copy`, other option is `breakpoints`
#' @param maxloherror Maximum value for LOH error rate
#' @param mincells Minimum cluster size used for phasing
#' @param overwritemincells Overwrite the the number of cells to use for clustering/phasing
#' @param cluster_per_chr Whether to cluster per chromosome to rephase alleles or not
#' @param filterhaplotypes filter out haplotypes present in less than X fraction, default is 0.1
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
#' sim_data <- simulate_data_cohort(
#'   clone_num = c(20, 20),
#'   clonal_events = list(
#'     list("1" = c(2, 0), "5" = c(3, 1)),
#'     list("2" = c(6, 3), "3" = c(1, 0))
#'   ),
#'   loherror = 0.02,
#'   coverage = 100
#' )
#'
#' results <- callHaplotypeSpecificCN(sim_data$CNbins, sim_data$haplotypes)
#' @md
#' @export
callHaplotypeSpecificCN <- function(CNbins,
                                    haplotypes,
                                    eps = 1e-12,
                                    loherror = 0.02,
                                    maxCN = NULL,
                                    selftransitionprob = 0.95,
                                    progressbar = TRUE,
                                    ncores = 1,
                                    phasebyarm = FALSE,
                                    minfrac = 0.7,
                                    likelihood = "binomial",
                                    minbins = 100,
                                    minbinschr = 10,
                                    phased_haplotypes = NULL,
                                    clustering_method = "copy",
                                    maxloherror = 0.035,
                                    mincells = 5,
                                    overwritemincells = NULL,
                                    cluster_per_chr = TRUE,
                                    viterbiver = "cpp", 
                                    filterhaplotypes = 0.1) {
  if (!clustering_method %in% c("copy", "breakpoints")) {
    stop("Clustering method must be one of copy or breakpoints")
  }
  
  if (!likelihood %in% c("binomial", "betabinomial", "auto")) {
    stop("Likelihood model for HMM emission model must be one of binomial, betabinomial or auto",
         call. = FALSE
    )
  }
  
  if (likelihood == "betabinomial" | likelihood == "auto") {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
           call. = FALSE
      )
    }
  }
  
  if (is.null(maxCN)) {
    maxCN <- max(CNbins$state)
  }
  
  nhaplotypes <- haplotypes %>% 
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(n = sum(totalcounts)) %>% 
    dplyr::pull(n) %>% 
    mean(.)
  
  haplotype_counts <- haplotypes %>% 
    dplyr::group_by(chr, start, end, hap_label) %>% 
    dplyr::summarize(totalcounts = sum(totalcounts), ncells = length(unique(cell_id))) %>% 
    dplyr::ungroup()
  
  haplotypes <- filter_haplotypes(haplotypes, filterhaplotypes)
  
  nhaplotypes_filt <- haplotypes %>% 
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(n = sum(totalcounts)) %>% 
    dplyr::pull(n) %>% 
    mean(.)
  
  cnbaf <- combineBAFCN(
    haplotypes = haplotypes,
    CNbins = CNbins,
    minbinschr = minbinschr,
    minbins = minbins
  )
  
  ascn <- .callHaplotypeSpecificCN_(cnbaf,
                                    eps = eps,
                                    loherror = loherror,
                                    maxCN = maxCN,
                                    selftransitionprob = selftransitionprob,
                                    progressbar = progressbar,
                                    ncores = ncores,
                                    viterbiver = viterbiver
  )
  
  infloherror <- ascn %>%
    dplyr::filter(state_phase == "A-Hom") %>%
    dplyr::summarise(err = weighted.mean(x = BAF, w = totalcounts, na.rm = TRUE)) %>%
    dplyr::pull(err)
  infloherror <- min(infloherror, maxloherror) # ensure loh error rate is < maxloherror
  
  if (likelihood == "betabinomial" | likelihood == "auto") {
    bbfit <- fitBB(ascn)
    if (bbfit$taronesZ > 5) {
      likelihood <- "betabinomial"
      message(paste0("Tarones Z-score: ", round(bbfit$taronesZ, 3), ", using ", likelihood, " model for inference."))
    } else {
      likelihood <- "binomial"
      message(paste0("Tarones Z-score: ", round(bbfit$taronesZ, 3), ", using ", likelihood, " model for inference."))
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
  
  ascn$balance <- ifelse(ascn$phase == "Balanced", 0, 1)
  
  if (is.null(phased_haplotypes)) {
    chrlist <- proportion_imbalance(ascn,
                                    haplotypes,
                                    phasebyarm = phasebyarm,
                                    minfrac = minfrac,
                                    mincells = mincells,
                                    clustering_method = clustering_method,
                                    overwritemincells = overwritemincells,
                                    cluster_per_chr = cluster_per_chr
    )
    propdf <- chrlist$propdf
    chrlist <- chrlist$chrlist
    phased_haplotypes <- phase_haplotypes_bychr(
      haplotypes = haplotypes,
      chrlist = chrlist,
      phasebyarm = phasebyarm
    )
  } else {
    chrlist <- NULL
    propdf <- NULL
  }
  
  cnbaf <- combineBAFCN(
    haplotypes = haplotypes,
    CNbins = CNbins,
    phased_haplotypes = phased_haplotypes,
    minbinschr = minbinschr,
    minbins = minbins
  )
  
  nhaplotypes_final <- cnbaf %>% 
    dplyr::group_by(cell_id) %>% 
    dplyr::summarize(n = sum(totalcounts)) %>% 
    dplyr::pull(n) %>% 
    mean(.)
  
  message(paste0("Total fraction of haplotype counts retained: ", round(nhaplotypes_final / nhaplotypes, 4), " (Total counts (Millions): ", round(nhaplotypes, 3) / 1e6, ", after filtering: ", round(nhaplotypes_filt, 3) / 1e6, ", after phasing:", round(nhaplotypes_final, 3) / 1e6, ")"))
  
  hscn_data <- .callHaplotypeSpecificCN_(cnbaf,
                                         eps = eps,
                                         loherror = infloherror,
                                         maxCN = maxCN,
                                         selftransitionprob = selftransitionprob,
                                         progressbar = progressbar,
                                         ncores = ncores,
                                         likelihood = likelihood,
                                         rho = bbfit$rho,
                                         viterbiver = viterbiver
  )
  
  # Output
  out <- list()
  class(out) <- "hscn"
  
  haplotype_counts <- dplyr::left_join(haplotype_counts, phased_haplotypes, by = c("chr", "start", "end", "hap_label")) %>% 
    dplyr::mutate(phase = ifelse(is.na(phase), "removed", phase)) %>% as.data.frame(.)
  
  out[["data"]] <- hscn_data %>% as.data.frame()
  out[["phasing"]] <- chrlist
  out[["haplotype_phasing"]] <- phased_haplotypes
  out[["haplotype_stats"]] <- haplotype_counts
  out[["haplotypesretained"]] <- propdf
  out[["loherror"]] <- infloherror
  out[["eps"]] <- eps
  out[["likelihood"]] <- bbfit
  out[["qc_summary"]] <- qc_summary(hscn_data)
  out[["qc_per_cell"]] <- qc_per_cell(hscn_data)
  
  out <- fix_assignments(out)
  
  return(out)
}

#' @export
fix_assignments <- function(hscn) {
  if (hscn$likelihood$likelihood == "binomial") {
    suppressWarnings(hscn_data <- hscn$data %>%
                       as.data.table() %>%
                       .[, rlid := data.table::rleid(state_AS_phased),
                         by = .(cell_id, chr)
                         ] %>%
                       .[, nbins := .N, by = c("rlid", "chr", "cell_id")] %>%
                       .[, Majprev := shift(Maj, fill = 0)] %>%
                       .[, Minprev := shift(Min, fill = 0)] %>%
                       .[, Majafter := shift(Maj, fill = 0, type = "lead")] %>%
                       .[, Minafter := shift(Min, fill = 0, type = "lead")] %>%
                       .[, pMinafter := ifelse(is.nan(Minafter / state),
                                               0.0,
                                               Minafter / state
                       )] %>% # nan check to stop 0/0
                       .[, pMinafter := fifelse(
                         pMinafter == 0.0,
                         pMinafter + hscn$loherror,
                         pMinafter
                       )] %>%
                       .[, pMinafter := fifelse(
                         pMinafter == 1.0,
                         pMinafter - hscn$loherror,
                         pMinafter
                       )] %>%
                       .[, pMinprev := ifelse(is.nan(Minprev / state),
                                              0.0, Minprev / state
                       )] %>% # nan check to stop 0/0
                       .[, pMinprev := fifelse(
                         pMinprev == 0.0,
                         pMinprev + hscn$loherror,
                         pMinprev
                       )] %>%
                       .[, pMinprev := fifelse(
                         pMinprev == 1.0,
                         pMinprev - hscn$loherror,
                         pMinprev
                       )] %>%
                       .[, LLafter := dbinom(alleleA, alleleA + alleleB, p = pMinafter)] %>%
                       .[, LLafter := fifelse(is.na(LLafter), 0, LLafter)] %>%
                       .[, LLprev := dbinom(alleleA, alleleA + alleleB, p = pMinprev, )] %>%
                       .[, LLprev := fifelse(is.na(LLprev), 0, LLafter)] %>%
                       .[, state_min := fifelse(
                         nbins == 1 & LLafter < LLprev,
                         Minprev,
                         Min
                       )] %>%
                       .[, state_min := fifelse(
                         nbins == 1 & LLafter >= LLprev,
                         Minafter,
                         Min
                       )] %>%
                       .[, Maj := state - state_min] %>%
                       .[, Min := state_min] %>%
                       .[, Min := fifelse(Min < 0, 0, Min)] %>%
                       .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
                       .[, Min := fifelse(Min > state, state, Min)] %>%
                       .[, Maj := fifelse(Maj > state, state, Maj)] %>%
                       add_states() %>%
                       dplyr::select(
                         -LLafter, -LLprev, -nbins, -pMinafter,
                         -pMinprev, -Minafter, -Minprev, -Majafter,
                         -Majprev, -rlid
                       ))
    
    hscn_data <- hscn_data %>%
      as.data.table() %>%
      .[, rlid := data.table::rleid(state_AS_phased), by = .(cell_id, chr)] %>%
      .[, alleleAtot := sum(alleleA), by = c("rlid", "chr", "cell_id")] %>%
      .[, alleleBtot := sum(alleleB), by = c("rlid", "chr", "cell_id")] %>%
      .[, pMin := ifelse(is.nan(Min / state),
                         0.0,
                         Min / state
      )] %>%
      # nan check to stop 0/0
      .[, pMin := fifelse(
        pMin == 0.0,
        pMin + hscn$loherror,
        pMin
      )] %>%
      .[, pMin := fifelse(
        pMin == 1.0,
        pMin - hscn$loherror,
        pMin
      )] %>%
      .[, pMaj := ifelse(is.nan(Maj / state),
                         0.0,
                         Maj / state
      )] %>%
      .[, pMaj := fifelse(
        pMaj == 0.0,
        pMaj + hscn$loherror,
        pMaj
      )] %>%
      .[, pMaj := fifelse(
        pMaj == 1.0,
        pMaj - hscn$loherror,
        pMaj
      )] %>%
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
      dplyr::select(
        -LLassigned, -LLother,
        -alleleBtot, -alleleAtot, -pMin, -pMaj, -rlid
      )
  } else {
    suppressWarnings(hscn_data <- hscn$data %>%
                       as.data.table() %>%
                       .[, rlid := data.table::rleid(state_AS_phased), by = .(cell_id, chr)] %>%
                       .[, nbins := .N, by = c("rlid", "chr", "cell_id")] %>%
                       .[, Majprev := shift(Maj, fill = 0)] %>%
                       .[, Minprev := shift(Min, fill = 0)] %>%
                       .[, Majafter := shift(Maj, fill = 0, type = "lead")] %>%
                       .[, Minafter := shift(Min, fill = 0, type = "lead")] %>%
                       .[, pMinafter := ifelse(is.nan(Minafter / state),
                                               0.0,
                                               Minafter / state
                       )] %>% # nan check to stop 0/0
                       .[, pMinafter := fifelse(
                         pMinafter == 0.0,
                         pMinafter + hscn$loherror,
                         pMinafter
                       )] %>%
                       .[, pMinafter := fifelse(
                         pMinafter == 1.0,
                         pMinafter - hscn$loherror,
                         pMinafter
                       )] %>%
                       .[, pMinprev := ifelse(is.nan(Minprev / state),
                                              0.0,
                                              Minprev / state
                       )] %>% # nan check to stop 0/0
                       .[, pMinprev := fifelse(
                         pMinprev == 0.0,
                         pMinprev + hscn$loherror,
                         pMinprev
                       )] %>%
                       .[, pMinprev := fifelse(
                         pMinprev == 1.0,
                         pMinprev - hscn$loherror,
                         pMinprev
                       )] %>%
                       .[, LLafter := VGAM::dbetabinom(alleleA, alleleA + alleleB, rho = hscn$likelihood$rho, p = pMinafter)] %>%
                       .[, LLafter := fifelse(
                         is.na(LLafter),
                         0,
                         LLafter
                       )] %>%
                       .[, LLprev := VGAM::dbetabinom(alleleA, alleleA + alleleB, rho = hscn$likelihood$rho, p = pMinprev)] %>%
                       .[, LLprev := fifelse(
                         is.na(LLprev),
                         0,
                         LLafter
                       )] %>%
                       .[, state_min := fifelse(
                         nbins == 1 & LLafter < LLprev,
                         Minprev,
                         Min
                       )] %>%
                       .[, state_min := fifelse(
                         nbins == 1 & LLafter >= LLprev,
                         Minafter,
                         Min
                       )] %>%
                       .[, Maj := state - state_min] %>%
                       .[, Min := state_min] %>%
                       .[, Min := fifelse(Min < 0, 0, Min)] %>%
                       .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
                       .[, Min := fifelse(Min > state, state, Min)] %>%
                       .[, Maj := fifelse(Maj > state, state, Maj)] %>%
                       add_states() %>%
                       dplyr::select(
                         -LLafter, -LLprev, -nbins, -pMinafter,
                         -pMinprev, -Minafter, -Minprev, -Majafter,
                         -Majprev, -rlid
                       ))
    
    hscn_data <- hscn_data %>%
      as.data.table() %>%
      .[, rlid := data.table::rleid(state_AS_phased)] %>%
      .[, alleleAtot := sum(alleleA), by = "rlid"] %>%
      .[, alleleBtot := sum(alleleB), by = "rlid"] %>%
      .[, pMin := ifelse(is.nan(Min / state), 0.0, Min / state)] %>%
      .[, pMin := fifelse(pMin == 0.0, pMin + hscn$loherror, pMin)] %>%
      .[, pMin := fifelse(pMin == 1.0, pMin - hscn$loherror, pMin)] %>%
      .[, pMaj := ifelse(is.nan(Maj / state), 0.0, Maj / state)] %>%
      .[, pMaj := fifelse(pMaj == 0.0, pMaj + hscn$loherror, pMaj)] %>%
      .[, pMaj := fifelse(pMaj == 1.0, pMaj - hscn$loherror, pMaj)] %>%
      .[, LLassigned := VGAM::dbetabinom(alleleAtot, alleleAtot + alleleBtot, rho = hscn$likelihood$rho, p = pMin)] %>%
      .[, LLother := VGAM::dbetabinom(alleleAtot, alleleAtot + alleleBtot, rho = hscn$likelihood$rho, p = pMaj)] %>%
      .[, state_min := fifelse(LLother < LLassigned, Maj, Min)] %>%
      .[, Maj := state - state_min] %>%
      .[, Min := state_min] %>%
      .[, Min := fifelse(Min < 0, 0, Min)] %>%
      .[, Maj := fifelse(Maj < 0, 0, Maj)] %>%
      .[, Min := fifelse(Min > state, state, Min)] %>%
      .[, Maj := fifelse(Maj > state, state, Maj)] %>%
      add_states() %>%
      dplyr::select(-LLassigned, -LLother, -alleleBtot, -alleleAtot, -pMin, -pMaj, -rlid)
  }
  
  hscn[["data"]] <- hscn_data %>% as.data.frame()
  hscn[["qc_summary"]] <- qc_summary(hscn_data)
  hscn[["qc_per_cell"]] <- qc_per_cell(hscn_data)
  return(hscn)
}

myphasedist <- function(A1, A2, B1, B2) {
  x <- sum(A1 != A2, na.rm = T) +
    sum(B1 != B2, na.rm = T)
  return(x)
}

#' @export
getphase <- function(A, B) {
  nbins <- dim(A)[1]
  
  Df <- data.frame(Dnon = 0, Dswa = 0)
  Bf <- data.frame(Bnon = -1, Bswa = -1)
  
  for (i in 2:nbins) {
    Dnonnon <- Df$Dnon[i - 1] + myphasedist(A[i, ], A[i - 1, ], B[i, ], B[i - 1, ])
    Dswanon <- Df$Dswa[i - 1] + myphasedist(A[i, ], B[i - 1, ], B[i, ], A[i - 1, ])
    Dnon <- min(Dnonnon, Dswanon)
    Bnon <- ifelse(Dnonnon <= Dswanon, 0, 1)
    
    Dnonswa <- Df$Dnon[i - 1] + myphasedist(B[i, ], A[i - 1, ], A[i, ], B[i - 1, ])
    Dswaswa <- Df$Dswa[i - 1] + myphasedist(B[i, ], B[i - 1, ], A[i, ], A[i - 1, ])
    Dswa <- min(Dnonswa, Dswaswa)
    Bswa <- ifelse(Dnonswa <= Dswaswa, 0, 1)
    
    Df <- dplyr::bind_rows(Df, data.frame(Dnon = Dnon, Dswa = Dswa))
    Bf <- dplyr::bind_rows(Bf, data.frame(Bnon = Bnon, Bswa = Bswa))
  }
  
  phasebin <- c()
  if (Df$Dnon[nbins] < Df$Dswa[nbins]) {
    phasebin <- c(phasebin, 0)
    prev <- Bf$Bnon[nbins]
  } else {
    phasebin <- c(phasebin, 1)
    prev <- Bf$Bswa[nbins]
  }
  
  for (i in (nbins - 1):1) {
    phasebintemp <- ifelse(prev == 0, 0, 1)
    phasebin <- c(phasebin, phasebintemp)
    prev <- Bf[, prev + 1][i]
  }
  
  phasebin <- !phasebin
  f <- table(phasebin)
  
  if (length(f) == 2) {
    if (f[2] > f[1]) {
      phasebin <- !phasebin
    }
  }
  
  phase <- unlist(lapply(phasebin, function(x) ifelse(x, "switch", "stick")))
  
  return(list(phase = rev(phase), Df = Df, Bf = Bf, phasebin = rev(phasebin)))
}

rephasehaplotypes <- function(phase_df, haplotype_phasing) {
  return(dplyr::left_join(haplotype_phasing, phase_df) %>%
           dplyr::mutate(phase = dplyr::case_when(
             phase == "allele0" & phasing == "stick" ~ "allele0",
             phase == "allele1" & phasing == "stick" ~ "allele1",
             phase == "allele0" & phasing == "switch" ~ "allele1",
             phase == "allele1" & phasing == "switch" ~ "allele0"
           )) %>%
           dplyr::select(chr, start, end, hap_label, phase))
}

phasing_minevents <- function(cndat, chromosomes) {
  phase_df <- data.frame()
  Amat <- createCNmatrix(cndat, field = "Maj")
  Bmat <- createCNmatrix(cndat, field = "Min")
  for (mychr in unique(cndat$chr)) {
    message(paste0("Phasing chromosome ", mychr))
    coords <- Amat %>%
      dplyr::filter(chr == mychr) %>%
      dplyr::select(chr, start, end)
    A <- Amat %>%
      dplyr::filter(chr == mychr) %>%
      dplyr::select(-chr, -start, -end, -width)
    B <- Bmat %>%
      dplyr::filter(chr == mychr) %>%
      dplyr::select(-chr, -start, -end, -width)
    if (mychr %in% chromosomes) {
      phase <- getphase(A, B)
      coords <- dplyr::bind_cols(coords, phase$Df)
      coords <- dplyr::bind_cols(coords, phase$Bf)
      coords$phasing <- phase$phase
    } else {
      coords$Dnon <- 1
      coords$Dswa <- 1
      coords$Bnon <- 1
      coords$Bswa <- 1
      coords$phasing <- "stick"
    }
    message(paste0("\tFraction of bins that will be switched: ", round(sum(coords$phasing == "switch") / dim(coords)[1], 2)))
    phase_df <- dplyr::bind_rows(phase_df, coords)
  }
  
  return(phase_df %>% dplyr::select(chr, start, end, phasing) %>% as.data.frame())
}

phasing_LOH <- function(cndat, chromosomes, cutoff = 0.9, ncells = 1) {
  phase_df <- data.frame()
  
  for (mychr in unique(cndat$chr)) {
    message(paste0("Phasing chromosome ", mychr))
    
    # find cells that have whole chromosome LOH
    cells <- cndat %>%
      as.data.table() %>%
      .[chr == mychr] %>%
      .[, list(
        LOH = sum(LOH == "LOH") / .N,
        # calculate distance between raw BAF and predicted state, we'll remove any cells that are far from this where the HMM may have failed
        mBAF = mean(BAF - (Min / state), na.rm = T)
      ),
      by = "cell_id"
      ] %>%
      .[order(LOH, decreasing = TRUE)] %>%
      .[LOH > 0.9 & abs(mBAF) < 0.05] %>%
      .$cell_id
    
    if (length(cells) > ncells) {
      BAFm <- cndat %>%
        as.data.table() %>%
        .[chr == mychr] %>%
        .[cell_id %in% cells] %>%
        .[, BAFm := ifelse(BAF > 0.5, 1 - BAF, BAF)] %>%
        .$BAFm %>%
        mean(.)
    } else {
      BAF <- 0.0
    }
    
    message(paste0("\tNumber of cells with whole chromosome LOH = ", length(cells)))
    
    if ((length(cells) > ncells) & (mychr %in% chromosomes)) {
      
      # create matrix of allele phases
      phasemat <- cndat %>%
        as.data.table() %>%
        .[cell_id %in% cells] %>%
        .[chr == mychr] %>%
        createCNmatrix(., field = "phase")
      
      coords <- phasemat %>% dplyr::select(chr, start, end)
      
      # find the average phase across these cells
      coords$averagephase <- apply(phasemat %>%
                                     dplyr::select(-chr, -start, -end, -width), 1, Mode)
      
      # find switches
      coords <- coords %>%
        dplyr::mutate(phasing = dplyr::case_when(
          averagephase == "A" ~ "stick",
          averagephase == "B" ~ "switch",
          averagephase == "Balanced" ~ "stick"
        ))
      message(paste0("\tFraction of bins that will be switched: ", round(sum(coords$phasing == "switch") / dim(coords)[1], 2)))
    } else {
      phasemat <- cndat %>%
        as.data.table() %>%
        .[chr == mychr] %>%
        createCNmatrix(., field = "phase")
      
      coords <- phasemat %>% dplyr::select(chr, start, end)
      coords$averagephase <- "A"
      coords <- coords %>%
        dplyr::mutate(phasing = dplyr::case_when(
          averagephase == "A" ~ "stick",
          averagephase == "B" ~ "switch",
          averagephase == "Balanced" ~ "stick"
        ))
    }
    
    phase_df <- dplyr::bind_rows(phase_df, coords)
  }
  
  return(phase_df %>% dplyr::select(chr, start, end, phasing) %>% as.data.frame(.))
}

#' Rephase haplotype copy number
#'
#' This function implements 2 rephasing algorithms. The first `mindist`, implementthe dynamic programming algorithm to rephase haplotype copy number first described in CHISEL.
#' The objective is to find the phase that minimizes the number of copy number events. The second `LOH` finds cells with whole chromosome losses and assumes this was a single
#' event and rephases all the bins relative to this.
#'
#' @param cn either a `hscn` object from `callHaplotypeSpecificCN` or a dataframe with haplotype specific copy number ie the `data` slot in an `hscn` object
#' @param chromosomes vector specifying which chromosomes to phase, default is NULL whereby all chromosomes are phased
#' @param method either `mindist` or `LOH`
#' @param ncells default 1
#'
#' @return Either a new `hscn` object or a dataframe with rephased bins depdending on the input
#'
#' @md
#' @export
rephasebins <- function(cn,
                        chromosomes = NULL,
                        method = "mindist",
                        whole_chr_cutoff = 0.9,
                        ncells = 1,
                        clusterfirst = FALSE,
                        cl = NULL) {
  if (!method %in% c("mindist", "LOH")) {
    stop(paste0("method must be one of mindist or LOH"))
  }
  
  
  if (is.hscn(cn)) {
    cndat <- cn$data
  } else {
    cndat <- cn
  }
  
  if (is.null(chromosomes)) {
    chromosomes <- unique(cndat$chr)
  }
  
  if (method == "mindist") {
    if (clusterfirst){
      if (is.null(cl)){
        ncells <- length(unique(cndat$cell_id))
        cl <- umap_clustering(cndat,
                              minPts = max(round(0.05 * ncells), 2),
                              field = "copy")
        cl <- cl$clustering
      }
      cndatclone <- consensuscopynumber(cndat, cl = cl)
      phase_df <- phasing_minevents(cndatclone, chromosomes)
    } else {
      phase_df <- phasing_minevents(cndat, chromosomes)
    }
  } else if (method == "LOH") {
    phase_df <- phasing_LOH(cndat,
                            chromosomes,
                            cutoff = whole_chr_cutoff,
                            ncells = ncells
    )
  }
  
  #switch the bins
  newhscn <- switchbins(cndat, phase_df)
  
  if (is.hscn(cn)) {
    cn$haplotype_phasing <-
      rephasehaplotypes(phase_df, cn$haplotype_phasing) %>% as.data.frame()
    cn$data <- newhscn %>% orderdf(.) %>%  as.data.frame()
    cn$qc_summary <- qc_summary(cn$data)
    x <- cn
  } else {
    x <- list(newhscn = newhscn %>% orderdf(.), phase = phase_df)
  }
  
  return(x)
}

#' @export
switchbins <- function(cn, phase_df){
  
  if (is.hscn(cn)) {
    cndat <- cn$data
  } else {
    cndat <- cn
  }
  
  newhscn <- as.data.table(cndat)[as.data.table(phase_df %>% dplyr::select(chr, start, end, phasing)),
                                  on = c("chr", "start", "end")
                                  ] %>%
    .[, Min := ifelse(phasing == "switch", Maj, Min)] %>%
    .[, Maj := ifelse(phasing == "switch", state - Min, Maj)] %>%
    .[, alleleB := ifelse(phasing == "switch", alleleA, alleleB)] %>%
    .[, alleleA := ifelse(phasing == "switch",
                          totalcounts - alleleB,
                          alleleA
    )] %>%
    .[, BAF := alleleB / totalcounts] %>%
    add_states()
  
  newhscn <- newhscn[, phasing := NULL]
  
  if (is.hscn(cn)) {
    cn$haplotype_phasing <-
      rephasehaplotypes(phase_df, cn$haplotype_phasing) %>% as.data.frame()
    cn$data <- newhscn %>% as.data.frame()
    cn$qc_summary <- qc_summary(cn$data)
    x <- cn
  } else {
    x <- newhscn
  }
  
  return(x)
}