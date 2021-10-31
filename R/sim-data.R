getstate <- function(base_ploidy) {
  x <- base_ploidy
  while (x == base_ploidy) {
    x <- round(rgamma(1, base_ploidy + 1, rate = 1))
  }
  if (x == 0) {
    x <- 1
  }
  return(x)
}

#' @export
simulate_cell <- function(nchr = 2,
                          hlamps = 2,
                          coverage = 10,
                          clonal_events = NULL,
                          base_ploidy = 2,
                          maxCN = 8,
                          copysd = 0.25,
                          likelihood = "binomial",
                          sampling_dist = "poisson",
                          rho = 0.0,
                          loherror = 0.01) {
  if (likelihood == "betabinomial") {
    if (!requireNamespace("VGAM", quietly = TRUE)) {
      stop("Package \"VGAM\" needed to use the beta-binomial model. Please install it.",
        call. = FALSE
      )
    }
  }

  data("dlpbins", envir = environment())
  bins <- dlpbins
  bins$state <- base_ploidy
  bins$state_min <- round(base_ploidy / 2)
  chrvec <- setdiff(unique(bins$chr), names(clonal_events))

  chromosomes <- list()
  if (nchr > 0) {
    for (i in 1:nchr) {
      chr_temp <- sample(chrvec, 1)
      totstate <- getstate(base_ploidy)
      Bstate <- round(runif(1, 0, floor(totstate / 2)))
      chromosomes[[chr_temp]] <- c(totstate, Bstate)
    }
  }

  chromosomes <- c(chromosomes, clonal_events)

  for (i in 1:length(chromosomes)) {
    bins <- bins %>%
      dplyr::mutate(state = ifelse(chr == names(chromosomes)[i], chromosomes[[i]][1], state))
    bins <- bins %>%
      dplyr::mutate(state_min = ifelse(chr == names(chromosomes)[i], chromosomes[[i]][2], state_min))
  }

  bins$copy <- bins$state + rnorm(length(bins$state), 0, copysd)
  bins$cell_id <- paste("sample", "lib", rawToChar(as.raw(sample(c(65:90, 97:122), 10, replace = T))), sep = "-")
  bins <- dplyr::select(bins, cell_id, chr, start, end, state, copy, state_min)


  ascn <- bins %>%
    data.table::as.data.table() %>%
    .[, Maj := state - state_min] %>%
    .[, Min := state_min] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := data.table::fifelse(state_min == 0, "LOH", "NO")] %>%
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
    .[, state_BAF := data.table::fifelse(is.nan(state_BAF), 0.5, state_BAF)]

  CNbins <- dplyr::select(ascn, cell_id, chr, start, end, state, copy)

  if (sampling_dist == "poisson"){
    haps <- ascn %>%
      .[, totalcounts := rpois(1, coverage * state), by = .(chr, start, end)]
  } else {
    haps <- ascn %>%
      .[, totalcounts := round(state * rgamma(1, shape = 1.76, rate = 0.09)), by = .(chr, start, end)]
  }

  if (likelihood == "binomial") {
    haps <- haps %>%
      .[, prob := fifelse(Min == 0, (Min / state) + loherror, Min / state)] %>%
      .[, alleleB := rbinom(n = 1, size = totalcounts, p = prob),
        by = .(chr, start, end)
      ] %>%
      .[, alleleA := totalcounts - alleleB] %>%
      .[, hap_label := 1:.N]
  } else {
    haps <- haps %>%
      .[, prob := fifelse(Min == 0, (Min / state) + loherror, Min / state)] %>%
      .[, alleleB := VGAM::rbetabinom(n = 1, size = totalcounts, p = prob, rho = rho), by = .(chr, start, end)] %>%
      .[, alleleA := totalcounts - alleleB] %>%
      .[, hap_label := 1:.N]
  }

  ascn$BAF <- haps$alleleB / haps$totalcounts

  haps <- haps %>%
    dplyr::select(cell_id, chr, start, end, hap_label, alleleA, alleleB, totalcounts)

  return(list(ascn = ascn, haplotypes = haps, CNbins = CNbins))
}

#' @export
simulate_cells <- function(ncells,
                           nchr = 2,
                           hlamps = 2,
                           coverage = 10,
                           clonal_events = NULL,
                           base_ploidy = 2,
                           maxCN = 8,
                           copysd = 0.25,
                           likelihood = "binomial",
                           sampling_dist = "poisson",
                           rho = 0.0,
                           loherror = 0.01) {
  ascn <- data.frame()
  CNbins <- data.frame()
  haps <- data.frame()
  for (i in 1:ncells) {
    print(i)
    cell <- simulate_cell(
      nchr = nchr,
      hlamps = hlamps,
      coverage = coverage,
      clonal_events = clonal_events,
      base_ploidy = base_ploidy,
      maxCN = maxCN,
      copysd = copysd,
      likelihood = likelihood,
      sampling_dist = sampling_dist,
      rho = rho,
      loherror = loherror
    )
    ascn <- dplyr::bind_rows(ascn, cell$ascn)
    CNbins <- dplyr::bind_rows(CNbins, cell$CNbins)
    haps <- dplyr::bind_rows(haps, cell$haplotypes)
  }

  switches <- dplyr::distinct(haps, chr, start, end, hap_label) %>%
    as.data.table() %>%
    .[, switch := sample(c(TRUE, FALSE), .N, TRUE), by = "hap_label"]

  haps <- haps %>%
    dplyr::left_join(switches) %>%
    as.data.table() %>%
    .[, allele1 := data.table::fifelse(switch == TRUE, alleleA, alleleB)] %>%
    .[, allele0 := data.table::fifelse(switch == TRUE, alleleB, alleleA)]

  haps <- haps %>%
    dplyr::select(cell_id, chr, start, end, hap_label, allele1, allele0, totalcounts, switch)

  return(list(
    ascn = ascn %>% as.data.frame(.),
    haplotypes = haps %>% as.data.frame(.),
    CNbins = CNbins %>% as.data.frame(.)
  ))
}

#' @export
simulate_data_cohort <- function(clone_num = c(10),
                                 coverage = 10,
                                 clonal_events = list(list("1" = c(2, 0))),
                                 base_ploidy = 2,
                                 nchr = 0,
                                 hlamps = 0,
                                 maxCN = 8,
                                 copysd = 0.25,
                                 likelihood = "binomial",
                                 sampling_dist = "poisson",
                                 rho = 0.0,
                                 loherror = 0.01) {
  ascn <- data.frame()
  CNbins <- data.frame()
  haps <- data.frame()
  for (clone in 1:length(clone_num)) {
    for (i in 1:clone_num[clone]) {
      cell <- simulate_cell(
        nchr = nchr,
        hlamps = hlamps,
        coverage = coverage,
        clonal_events = clonal_events[[clone]],
        base_ploidy = base_ploidy,
        maxCN = maxCN,
        copysd = copysd,
        likelihood = likelihood,
        sampling_dist = sampling_dist,
        rho = rho,
        loherror = loherror
      )
      ascn <- dplyr::bind_rows(ascn, cell$ascn)
      CNbins <- dplyr::bind_rows(CNbins, cell$CNbins)
      haps <- dplyr::bind_rows(haps, cell$haplotypes)
    }
  }

  switches <- dplyr::distinct(haps, chr, start, end, hap_label) %>%
    as.data.table() %>%
    .[, switch := sample(c(TRUE, FALSE), .N, TRUE), by = "hap_label"]

  haps <- haps %>%
    dplyr::left_join(switches) %>%
    as.data.table() %>%
    .[, allele1 := data.table::fifelse(switch == TRUE, alleleA, alleleB)] %>%
    .[, allele0 := data.table::fifelse(switch == TRUE, alleleB, alleleA)]

  haps <- haps %>%
    dplyr::select(cell_id, chr, start, end, hap_label, allele1, allele0, totalcounts, switch)

  return(list(
    ascn = ascn %>% as.data.frame(.),
    haplotypes = haps %>% as.data.frame(.),
    CNbins = CNbins %>% as.data.frame(.)
  ))
}
