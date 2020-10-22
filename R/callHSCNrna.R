get_states_dna <- function(hscn, minf = 0.1, arms = NULL){

  if (is.null(arms)){
    possible_states <- hscn %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(state_AS_phased = Mode(state_AS_phased),
                       Maj = Mode(Maj),
                       Min = Mode(Min)) %>%
      dplyr::group_by(chr, arm, chrarm, state_AS_phased, Min, Maj) %>%
      dplyr::summarise(n = n(), f = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(chr, arm, chrarm) %>%
      dplyr::mutate(f = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(f > minf)
  } else {
    possible_states <- hscn %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = FALSE)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(arm = ifelse(chrarm %in% arms, arm, "")) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(state_AS_phased = Mode(state_AS_phased),
                       Maj = Mode(Maj),
                       Min = Mode(Min)) %>%
      dplyr::group_by(chr, arm, chrarm, state_AS_phased, Min, Maj) %>%
      dplyr::summarise(n = n(), f = sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(chr, arm, chrarm) %>%
      dplyr::mutate(f = n / sum(n)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(f > minf)
  }

  return(possible_states)
}

#' @export
assign_states <- function(haps,
                          hscn,
                          minf = 0.1,
                          shrinkage = FALSE,
                          loherror = 0.03,
                          mergelowcounts = TRUE,
                          arms = NULL){

  perchrlist <- per_arm_baf_mat(haps, mergelowcounts = mergelowcounts, arms = arms)
  bafperchr <- perchrlist$bafperchr
  possible_states <- get_states_dna(hscn, minf = minf, arms = unique(bafperchr$chrarm))

  data("hg19chrom_coordinates", envir=environment())

  if (shrinkage == FALSE){

    perchr <- dplyr::left_join(bafperchr, possible_states) %>%
      dplyr::arrange(cell_id, chrarm, state_AS_phased) %>%
      as.data.table() %>%
      .[, prob := Min / (Min + Maj)] %>%
      .[, prob := fifelse(prob == 0.0, prob + loherror, prob)] %>%
      .[, prob := fifelse(prob == 1.0, prob - loherror, prob)] %>%
      .[, L := dbinom(size = total, x = alleleB, prob = prob)] %>%
      .[, state := Min + Maj]
    perchr <- perchr[perchr[, .I[which.max(L)], by=.(chrarm, cell_id)]$V1] %>%
      add_states() %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      as.data.table() %>%
      .[, copy := state]

    perchr <- as.data.frame(perchr) %>%
      dplyr::select(-L, -n, -f, -prob) %>%
      as.data.frame()
  } else {

    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
           call. = FALSE)
    }

    # negative log likelihood of data given alpha; beta
    llfunc <- function(data){
      ll <- function(alpha, beta) {
        -sum(VGAM::dbetabinom.ab(data$alleleB, data$total, alpha, beta, log = TRUE))
      }
    }

    # identify MLE for alpha and beta per chromosome arm
    m <- lapply(X = unique(bafperchr$chrarm),
                function(x) {
                  m <- stats4::mle(llfunc(dplyr::filter(bafperchr, chrarm == x)),
                          start = list(alpha = 15, beta = 15), method = "Nelder-Mead")
                  m_coef <- stats4::coef(m)
                  return(data.frame(alpha = m_coef[["alpha"]],
                                    beta = m_coef[["beta"]],
                                    chrarm = x))
                  })

    # add mean of distribution
    m <- dplyr::bind_rows(m) %>%
      dplyr::mutate(mle_f = alpha / (alpha + beta))

    #use empirical bayes to estimate expected BAF
    bafperchr <- dplyr::left_join(bafperchr, m) %>%
      dplyr::mutate(newBAF = (alleleB + alpha) / (total + alpha + beta))

    perchr <- dplyr::left_join(bafperchr, possible_states) %>%
      dplyr::arrange(cell_id, chrarm, state_AS_phased) %>%
      as.data.table() %>%
      .[, prob := Min / (Min + Maj)] %>%
      .[, prob := fifelse(prob == 0.0, prob + loherror, prob)] %>%
      .[, prob := fifelse(prob == 1.0, prob - loherror, prob)] %>%
      .[, L := dbinom(size = round(total + alpha + beta), #empirical bayes correction
                      x = round(alleleB + alpha),
                      prob = prob)] %>%
      .[, state := Min + Maj]
    perchr <- perchr[perchr[, .I[which.max(L)], by=.(chrarm, cell_id)]$V1] %>%
      add_states() %>%
      .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)] %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      as.data.table() %>%
      .[, copy := state]

    perchr <- as.data.frame(perchr) %>%
      dplyr::select(-L, -n, -f, -prob) %>%
      as.data.frame()
  }

  return(perchr)
}

possible_states_df <- function(bafperchr, step = 0.25){

  possible_states <- expand.grid(chrarm = unique(bafperchr$chrarm),
              prob = seq(0.0, 1.0, step)) %>%
    dplyr::mutate(state = 1/step, Min = prob * state) %>%
    dplyr::mutate(Maj = state - Min) %>%
    dplyr::mutate(state_AS_phased = paste0(Maj, "|", Min))

  return(possible_states)
}

#' @export
assign_states_noprior <- function(haps,
                                  mergelowcounts = TRUE,
                                  shrinkage = FALSE,
                                  step = 0.25,
                                  loherror = 0.03,
                                  arms = NULL){

  perchrlist <- per_arm_baf_mat(haps, mergelowcounts = mergelowcounts, arms = arms)
  bafperchr <- perchrlist$bafperchr
  possible_states <- possible_states_df(bafperchr, step = step)

  data("hg19chrom_coordinates", envir=environment())

  if (shrinkage == FALSE){

    perchr <- dplyr::left_join(bafperchr, possible_states) %>%
      dplyr::arrange(cell_id, chrarm, state_AS_phased) %>%
      as.data.table() %>%
      .[, prob := Min / (Min + Maj)] %>%
      .[, prob := fifelse(prob == 0.0, prob + loherror, prob)] %>%
      .[, prob := fifelse(prob == 1.0, prob - loherror, prob)] %>%
      .[, L := dbinom(size = total, x = alleleB, prob = prob)] %>%
      .[, state := Min + Maj]

    perchr <- perchr[perchr[, .I[which.max(L)], by=.(chrarm, cell_id)]$V1] %>%
      add_states() %>%
      .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)] %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      as.data.table() %>%
      .[, copy := state]

    perchr <- as.data.frame(perchr) %>%
      #dplyr::select(-L, -prob) %>%
      as.data.frame()
  } else {

    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
           call. = FALSE)
    }

    # negative log likelihood of data given alpha; beta
    llfunc <- function(data){
      ll <- function(alpha, beta) {
        -sum(VGAM::dbetabinom.ab(data$alleleB, data$total, alpha, beta, log = TRUE))
      }
    }

    # identify MLE for alpha and beta per chromosome arm
    m <- lapply(X = unique(bafperchr$chrarm),
                function(x) {
                  m <- stats4::mle(llfunc(dplyr::filter(bafperchr, chrarm == x)),
                                   start = list(alpha = 15, beta = 15), method = "Nelder-Mead")
                  m_coef <- stats4::coef(m)
                  return(data.frame(alpha = m_coef[["alpha"]],
                                    beta = m_coef[["beta"]],
                                    chrarm = x))
                })

    # add mean of distribution
    m <- dplyr::bind_rows(m) %>%
      dplyr::mutate(mle_f = alpha / (alpha + beta))

    #use empirical bayes to estimate expected BAF
    bafperchr <- dplyr::left_join(bafperchr, m) %>%
      dplyr::mutate(newBAF = (alleleB + alpha) / (total + alpha + beta))

    perchr <- dplyr::left_join(bafperchr, possible_states) %>%
      dplyr::arrange(cell_id, chrarm, state_AS_phased) %>%
      as.data.table() %>%
      .[, prob := Min / (Min + Maj)] %>%
      .[, prob := fifelse(prob == 0.0, prob + loherror, prob)] %>%
      .[, prob := fifelse(prob == 1.0, prob - loherror, prob)] %>%
      .[, L := dbinom(size = round(total + alpha + beta), #empirical bayes correction
                      x = round(alleleB + alpha),
                      prob = prob)] %>%
      .[, state := Min + Maj]
    perchr <- perchr[perchr[, .I[which.max(L)], by=.(chrarm, cell_id)]$V1] %>%
      add_states() %>%
      .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)] %>%
      dplyr::left_join(hg19chrom_coordinates) %>%
      as.data.table() %>%
      .[, copy := state]

    perchr <- as.data.frame(perchr) %>%
      #dplyr::select(-L, -prob) %>%
      as.data.frame()
  }

  return(perchr)
}

#' @export
assign_states_dp <- function(bafperchr,
                             samples = 10,
                             alpha_0 = 1,
                             K = 20,
                             most_variable_chr = TRUE,
                             top_nchr = 5,
                             overwrite_chr = NULL,
                             removechr = c("chrX"),
                             filtercounts = 0){

  data("hg19chrom_coordinates", envir=environment())

  if (!requireNamespace("VIBER", quietly = TRUE)) {
    stop("Package \"VIBER\" needed to use the beta-binomial model.",
         call. = FALSE)
  }

  message('Generating count matrices...')
  baf_total <- bafperchr %>%
    dplyr::select(cell_id, chrarm, total) %>%
    tidyr::pivot_wider(names_from = "chrarm",
                       values_from = "total",
                       values_fill = list(total = 0),
                       names_prefix = "chr") %>%
    as.data.frame()
  row.names(baf_total) <- baf_total$cell_id
  baf_total = subset(baf_total, select = -c(cell_id))

  # make chr ~ cell_id matrix for B allele counts
  baf_counts <- bafperchr %>%
    dplyr::select(cell_id, chrarm, alleleB) %>%
    tidyr::pivot_wider(names_from = "chrarm",
                       values_from = "alleleB",
                       values_fill = list(alleleB = 0),
                       names_prefix = "chr") %>%
    as.data.frame()
  row.names(baf_counts) <- baf_counts$cell_id
  baf_counts = subset(baf_counts, select = -c(cell_id))

  if (filtercounts > 0){
    filtered_cells <- rowSums(baf_total) > filtercounts
    filtered_cells <- names(filtered_cells[filtered_cells == TRUE])
    baf_counts <- baf_counts[filtered_cells, ]
    baf_total <- baf_total[filtered_cells, ]
  }

  if (most_variable_chr){
    #keepchrs <- sort(sapply(baf_counts / baf_total, function(x) var(x, na.rm = TRUE)), decreasing = TRUE)[1:top_nchr]
    chr_names <- names(baf_counts)
    baf_counts_temp <- baf_counts[, setdiff(chr_names, removechr)]
    keepchrs <- sort(sapply(names(baf_counts_temp),
                              function(x) matrixStats::weightedVar(baf_counts[,x] / baf_total[,x],
                                                       w = baf_total[,x], na.rm = TRUE)),
                     decreasing = TRUE)[1:top_nchr]
    message(paste0("Top ", top_nchr, " most variable chromosomes are: ", paste0(names(keepchrs), collapse = ", ")))
    message("Using these chromosomes for clustering")
    keepchrs <- names(keepchrs)
  } else if (!is.null(overwrite_chr)) {
    keepchrs <- overwrite_chr
  } else{
    keepchrs <- paste0("chr", unique(bafperchr$chrarm))
  }

  baf_counts <- baf_counts[,keepchrs]
  baf_total <- baf_total[,keepchrs]
  print(dim(baf_counts))

  message('Fitting mixture model using VIBER...')
  fit = VIBER::variational_fit(
    baf_counts,
    baf_total,
    K = K,
    alpha_0 = alpha_0,
    samples = samples,
    q_init = "prior"
  )

  message('Filtering mixture components...')
  fit_filt <- VIBER::choose_clusters(fit,
                              binomial_cutoff = 0,
                              dimensions_cutoff = 0,
                              pi_cutoff = 0.01)

  message('Extract cluster means from VIBER object...')
  theta <- as.data.frame(fit_filt$theta_k)
  theta$chrarm <- row.names(theta)
  row.names(theta) <- NULL
  theta <- theta %>%
    tidyr::pivot_longer(-chrarm, names_to = "clone_id", values_to = "theta") %>%
    dplyr::mutate(chrarm = str_remove_all(chrarm, "chr"))

  message('Generating dataframe mapping cell_id to clone_id...')
  x <- data.frame(clone_id = fit_filt$labels$cluster.Binomial,
                  cell_id = row.names(baf_counts))

  message('Assign states to clones and chromosomes...')
  states <- bafperchr %>%
    dplyr::left_join(x, by = "cell_id") %>%
    dplyr::group_by(chrarm, clone_id) %>%
    dplyr::summarise(BAF = median(BAF),
              alleleA = sum(alleleA),
              alleleB = sum(alleleB),
              total = sum(total)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(theta, by = c("clone_id", "chrarm")) %>%
    dplyr::mutate(rounded = round(BAF / 0.25) * 0.25) %>%
    dplyr::mutate(state_phase = dplyr::case_when(
      rounded == 0.0 ~ "A-LOH",
      rounded == 0.25 ~ "A-Gained",
      rounded == 0.5 ~ "Balanced",
      rounded == 0.75 ~ "B-Gained",
      rounded == 1.0 ~ "B-LOH"
    )) %>%
    dplyr::select(chrarm, clone_id, state_phase)

  message("Format final dataframe...")
  bafperchr_new <- bafperchr %>%
    dplyr::select(chr, arm, chrarm, cell_id, alleleA, alleleB, total, BAF) %>%
    dplyr::left_join(x, by = "cell_id") %>%
    dplyr::left_join(states, by = c("chrarm", "clone_id")) %>%
    dplyr::select(-clone_id) %>%
    dplyr::left_join(hg19chrom_coordinates)


  return(list(viber_fit = fit_filt, clusters = x, ascn = bafperchr_new, usedchrs = keepchrs))
}

