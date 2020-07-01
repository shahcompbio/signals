plottinglist <- function(CNbins, xaxis_order = "genome_position", maxCN = 20){
  #arrange segments in order, generate segment index and reorder CN state factor

  if (!xaxis_order %in% c("bin", "genome_position")){
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }

  if (xaxis_order == "bin"){
  chridx <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))
  CNbins <- CNbins %>%
    dplyr::filter(!is.na(copy)) %>%
    dplyr::left_join(., chridx, by = c("chr")) %>%
    dplyr::arrange(cell_id, idx, start) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(copy = ifelse(copy > maxCN, maxCN, copy)) %>%
    dplyr::mutate(state = ifelse(state > maxCN, maxCN, state)) %>%
    dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>%
    dplyr::mutate(CNs = forcats::fct_reorder(ifelse(is.na(state), NA,
                                                    paste0("CN", state)), state))

  #get breaks - first index of each chromosome
  chrbreaks <- CNbins %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(idx)

  chrlabels <- gtools::mixedsort(unique(CNbins$chr))

  minidx = min(CNbins$idx)
  maxidx <- max(CNbins$idx)
  } else {
    binsize <- CNbins$end[1] - CNbins$start[1] + 1
    bins <- getBins(binsize = binsize) %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::mutate(idx = 1:dplyr::n())


    CNbins <- dplyr::full_join(bins, CNbins) %>%
      dplyr::filter(!is.na(copy)) %>%
      dplyr::filter(!is.na(state)) %>%
      dplyr::mutate(copy = ifelse(copy > maxCN, maxCN, copy)) %>%
      dplyr::mutate(state = ifelse(state > maxCN, maxCN, state)) %>%
      dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>%
      dplyr::mutate(CNs = forcats::fct_reorder(ifelse(is.na(state), NA,
                                                      paste0("CN", state)), state))

    #get breaks - first index of each chromosome
    chrbreaks <- bins %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::group_by(chr) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(idx)

    chrlabels <- gtools::mixedsort(unique(CNbins$chr))
    minidx <- min(bins$idx)
    maxidx <- max(bins$idx)
  }

  return(list(CNbins = CNbins, chrbreaks = chrbreaks, chrlabels = chrlabels, minidx = minidx, maxidx = maxidx))
}

#' @export
plot_umap <- function(clustering, bycol = NULL, alphavalue = 0.5, raster = FALSE){

  if (raster == TRUE){
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    g <- ggplot2::ggplot(clustering, ggplot2::aes(x = umap1, y = umap2)) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = bycol), alpha = alphavalue) +
      ggplot2::xlab("UMAP 1") +
      ggplot2::ylab("UMAP 2") +
      ggplot2::theme_bw() +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
  } else {
    g <- ggplot2::ggplot(clustering, ggplot2::aes(x = umap1, y = umap2)) +
      ggplot2::geom_point(ggplot2::aes_string(col = bycol), alpha = alphavalue) +
      ggplot2::xlab("UMAP 1") +
      ggplot2::ylab("UMAP 2") +
      ggplot2::theme_bw() +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 3))
  }
  return(g)
}

#' @export
squashy_trans <- function() {
  scales::trans_new(
    "squashy",
    function(x) tanh( 0.075 * x),
    function(x) atanh(x) / 0.075
  )
}

#' @export
plotCNprofile <- function(CNbins,
                         cellid = NULL,
                         chrfilt = NULL,
                         pointsize = 1,
                         alphaval = 0.9,
                         maxCN = 10,
                         cellidx = 1,
                         statecol = "state",
                         returnlist = FALSE,
                         raster = FALSE,
                         y_axis_trans = "identity",
                         xaxis_order = "genome_position"){

  if (!xaxis_order %in% c("bin", "genome_position")){
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }

  if (is.null(cellid)){
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  if (y_axis_trans == "squashy"){
    maxCN <- min(c(20, maxCN))
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)){
    message(paste0("Filtering for chromosome: ", chrfilt))
    CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN)

  if (raster == TRUE){
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
      ggplot2::scale_color_manual(name = "Copy number",
                                  breaks = names(statecolpal),
                                  labels = names(statecolpal),
                                  values = statecolpal,
                                  drop = FALSE) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN),trans = y_axis_trans) +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                     override.aes = list(alpha=1, size = 3, shape = 15))) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
  } else {
    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
      ggplot2::scale_color_manual(name = "Allele Specific CN",
                                  breaks = names(statecolpal),
                                  labels = names(statecolpal),
                                  values = statecolpal,
                                  drop = FALSE) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN), trans = y_axis_trans) +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                     override.aes = list(alpha=1, size = 3, shape = 15))) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
  }

  if (returnlist == TRUE){
    gCN <- list(CN = gCN, plist = pl)
  }

  return(gCN)
}


#' @export
plotCNprofileBAF <- function(cn,
                          cellid = NULL,
                          chrfilt = NULL,
                          pointsize = 1,
                          alphaval = 0.9,
                          maxCN = 10,
                          cellidx = 1,
                          BAFcol = "state_phase",
                          statecol = "state",
                          returnlist = FALSE,
                          raster = FALSE,
                          y_axis_trans = "identity",
                          xaxis_order = "genome_position"){

  if (!xaxis_order %in% c("bin", "genome_position")){
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }

  if (is.hscn(cn) | is.ascn(cn)){
    CNbins <- cn$data
  } else{
    CNbins <- cn
  }

  if (y_axis_trans == "squashy"){
    maxCN <- min(c(20, maxCN))
  }

  if (is.null(cellid)){
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  if (!"BAF" %in% names(CNbins)){
    stop("No BAF column in dataframe, first calculate the BAF per bin using combineBAFCN and then callAlleleSpecificCN")
  }

  if (BAFcol == "state_min"){
    BAFcolpal <- scCNminorallele_cols()
  }

  if (BAFcol == "state_phase"){
    BAFcolpal <- scCNphase_cols()
  }

  if (BAFcol == "state_AS"){
    BAFcolpal <- scCNAS_cols()
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)){
    message(paste0("Filtering for chromosome: ", chrfilt))
    CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN)

  if (raster == TRUE){
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
    gBAF <- pl$CNbins %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = BAF)) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = BAFcol), size = pointsize, alpha = alphaval) +
      ggplot2::scale_color_manual(name = "CN",
                                  breaks = names(BAFcolpal),
                                  labels = names(BAFcolpal),
                                  values = BAFcolpal) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.0)) +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("BAF") +
      ggplot2::ggtitle(cellid) +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      ggplot2::geom_hline(yintercept = 0.5, lty = 2, alpha = 0.5) +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank()) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha=1, size = 3, shape = 15)))

    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
      ggplot2::scale_color_manual(name = "Allele Specific CN",
                                  breaks = names(statecolpal),
                                  labels = names(statecolpal),
                                  values = statecolpal,
                                  drop = FALSE) +
      ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     legend.position = "none") +
      ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN), trans = y_axis_trans) +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot() +
      cowplot::background_grid(major = "x") +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                     override.aes = list(alpha=1, size = 3, shape = 15))) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
    } else {
      gBAF <- pl$CNbins %>%
        dplyr::mutate(state_min = paste0(state_min)) %>%
        ggplot2::ggplot(ggplot2::aes(x = idx, y = BAF)) +
        ggplot2::geom_point(ggplot2::aes_string(col = BAFcol), size = pointsize, alpha = alphaval) +
        ggplot2::scale_color_manual(name = "CN",
                                    breaks = names(BAFcolpal),
                                    labels = names(BAFcolpal),
                                    values = BAFcolpal) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       legend.position = "none") +
        ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.0)) +
        ggplot2::xlab("Chromosome") +
        ggplot2::ylab("BAF") +
        ggplot2::ggtitle(cellid) +
        cowplot::theme_cowplot() +
        cowplot::background_grid(major = "x") +
        ggplot2::geom_hline(yintercept = 0.5, lty = 2, alpha = 0.5) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank()) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom") +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha=1, size = 3, shape = 15)))

      gCN <- pl$CNbins %>%
        dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
        dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
        dplyr::mutate(state_min = paste0(state_min)) %>%
        ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
        ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
        ggplot2::scale_color_manual(name = "Allele Specific CN",
                                    breaks = names(statecolpal),
                                    labels = names(statecolpal),
                                    values = statecolpal,
                                    drop = FALSE) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       axis.ticks.y = ggplot2::element_blank(),
                       legend.position = "none") +
        ggplot2::scale_x_continuous(breaks = pl$chrbreaks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) + #,guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN), trans = y_axis_trans) +
        ggplot2::xlab("Chromosome") +
        ggplot2::ylab("Copy Number") +
        cowplot::theme_cowplot() +
        cowplot::background_grid(major = "x") +
        ggplot2::guides(colour = ggplot2::guide_legend(ncol = 6, byrow = TRUE,
                                                       override.aes = list(alpha=1, size = 3, shape = 15))) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")
    }

  g <- cowplot::plot_grid(gBAF, gCN, align = "v", ncol = 1, rel_heights = c(1, 1.2))

  if (returnlist == TRUE){
    p <- list(CN = gCN, BAF = gBAF, both = g, plist = pl)
  } else {
    p <- g
  }

  return(p)
}

#' @export
plotCNBAF <- function(cn, nfilt = 10^5, plottitle = "5Mb", pointsize = 0.1){
  if (is.hscn(cn) | is.ascn(cn)){
    CNbins <- cn$data
  } else{
    CNbins <- cn
  }
  maj <- seq(1, 11, 1)
  min <- seq(0, 11, 1)
  ASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS = paste0(state-min, "|", min)) %>%
    dplyr::mutate(Maj = state-min, Min = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(Maj >= 0, Min >= 0)

  CNbins <- CNbins %>%
    dplyr::filter(state < 10, copy < 10)

  if (nfilt < length(CNbins$chr)){
    CNbins <- CNbins %>%
      dplyr::sample_n(nfilt)
  }

  plottitle <- paste0((CNbins$end[1] - CNbins$start[1] + 1) / 1e6, "Mb bins")

  g <- CNbins %>%
    ggplot2::ggplot(ggplot2::aes(y = BAF, x = copy, col = paste0("CN", state))) +
    ggplot2::geom_point(size = pointsize, alpha = 0.2) +
    ggplot2::xlab("Corrected read counts") +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle(plottitle) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size=5))) +
    ggplot2::geom_text(data = ASstates, ggplot2::aes(x = state, y = cBAF, label = state_AS), col = "black") +
    ggplot2::scale_color_manual(name = "Copy number \n state",
                       breaks = paste0("CN", seq(0, max(CNbins$state, na.rm = TRUE), 1)),
                       labels = seq(0, max(CNbins$state, na.rm = TRUE), 1),
                       values = scCN_cols(paste0("CN", seq(0, max(CNbins$state, na.rm = TRUE), 1)))) +
    ggplot2::scale_x_continuous(breaks = seq(0, 10, 1), labels = seq(0, 10, 1), limits = c(0, 10))
  return(g)
}

#' @export
plotBAFperstate <- function(cn, minpts = 250, minfrac = 0.01, maxstate = 10, dens_adjust = 2.0){

  if (is.hscn(cn) | is.ascn(cn)){
    alleleCN <- cn$data
  } else{
    alleleCN <- cn
  }

  maj <- seq(0, max(alleleCN$state), 1)
  min <- seq(0, max(alleleCN$state), 1)

  allASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS_phased = paste0(state-min, "|", min)) %>%
    dplyr::mutate(Maj = state-min, Min = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(Maj >= 0, Min >= 0) %>%
    dplyr::filter(state <= maxstate)
  allASstates$cBAF[is.nan(allASstates$cBAF)] <- 0.0

  forplot <- alleleCN %>%
    dplyr::group_by(state_AS_phased) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(f = n / dplyr::n()) %>%
    dplyr::filter(state <= maxstate, n > minpts, f > minfrac)
  allASstates <- allASstates %>%
    dplyr::filter(state_AS_phased %in% unique(forplot$state_AS_phased))

  text_fraction <- dplyr::distinct(forplot, state_AS_phased, f) %>%
    dplyr::mutate(pct = paste0(100 * round(f, 2), "%"), y = 0.95)

  g <- forplot %>%
    dplyr::mutate(cncol = paste0("CN", state)) %>%
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(state_AS_phased, state),
                                 y = BAF)) +
    ggplot2::geom_violin(scale = "width", col = "white", ggplot2::aes(fill = cncol), adjust = dens_adjust) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, col = "white", ggplot2::aes(fill = cncol)) +
    ggplot2::scale_fill_manual(name = "Copy number \n state",
                                breaks = paste0("CN", seq(0, max(alleleCN$state, na.rm = TRUE), 1)),
                                labels = seq(0, max(alleleCN$state, na.rm = TRUE), 1),
                                values = scCN_cols(paste0("CN", seq(0, max(alleleCN$state, na.rm = TRUE), 1)))) +
    cowplot::theme_cowplot() +
    ggplot2::geom_crossbar(data = allASstates, ggplot2::aes(y = cBAF, ymin = cBAF, ymax = cBAF),
                           alpha = 0.2, size = 0.2) +
    ggplot2::geom_text(data = text_fraction, ggplot2::aes(x = state_AS_phased, y = y, label = pct)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(g)
}


plot_density_histogram <- function(dat, mystate, rho, nbins = 30, frac = "NA", loherror = 0.02){
  dat <- dat %>%
    dplyr::filter(state_AS_phased == mystate, totalcounts > 9)

  expBAF <- strsplit(mystate, "\\|")[[1]]
  expBAF <- as.numeric(expBAF[2]) / (as.numeric(expBAF[1]) + as.numeric(expBAF[2]))
  if (expBAF > 0.5) {
    expBAF <- expBAF - loherror
  } else if (expBAF < 0.5) {
    expBAF <- expBAF + loherror
  } else{
    expBAF <- expBAF
  }

  BAF_bb <- VGAM::rbetabinom(length(dat$totalcounts), size = dat$totalcounts,
                             prob = expBAF,
                             rho = rho) / dat$totalcounts
  dffit_bb <- data.frame(BAF = BAF_bb,
                         type = paste("Beta Binomial (rho = ", round(rho, 4), ")", sep=""))

  BAF_b <- VGAM::rbetabinom(length(dat$totalcounts), size = dat$totalcounts,
                            prob = expBAF,
                            rho = 0.0) / dat$totalcounts
  dffit_b <- data.frame(BAF = BAF_b,
                        type = paste("Binomial", sep=""))

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = BAF)) +
    ggplot2::geom_histogram(ggplot2::aes(y=(..density..)), bins = nbins, fill = "azure4",alpha = 0.8) +
    ggplot2::geom_line(data = dffit_bb, stat="density",ggplot2::aes(BAF, col = type), size = 0.5, adjust = 5) +
    ggplot2::geom_line(data = dffit_b, stat="density",ggplot2::aes(BAF, col = type), size = 0.5, adjust = 5, linetype = 2) +
    ggplot2::scale_colour_manual(values=c(ggplot2::alpha("deepskyblue4",0.6), ggplot2::alpha("firebrick4",0.6))) +
    ggplot2::theme_bw() +
    ggplot2::xlab("BAF") +
    ggplot2::ylab("Density") +
    ggplot2:: xlim(c(0.0, 1.0)) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none") +
    ggplot2:: ggtitle(paste0(mystate, " (", round(frac, 3) * 100, "%)"))

  return(g)
}

#' @export
plotBBfit <- function(hscn, nbins = 30, minfrac = 0.01){

  x <- table(hscn$data$state_AS_phased)
  x <- x / sum(x)
  x <- x[which(x > minfrac)]

  mydat <- dplyr::filter(hscn$data, Maj != 0, Min != 0, totalcounts > 9)

  x <- x[names(x) %in% unique(mydat$state_AS_phased)]

  gplots <- list()
  j <- 1
  for (i in names(x)){
    gplots[[j]] <- plot_density_histogram(hscn$data,
                                          mystate = i,
                                          rho = hscn$likelihood$rho,
                                          nbins = nbins,
                                          frac = x[[i]],
                                          loherror = hscn$loherror)
    j <- j + 1
  }

  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    gplots[[j-1]] + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12)) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) + ggplot2::theme(legend.position = 'right')
  )

  gplots[[j]] <- legend

  cowplot::plot_grid(plotlist = gplots, ncol = 3)
}

#' @export
plot_variance_state <- function(hscn, by_allele_specific_state = FALSE){
  if (by_allele_specific_state == TRUE){
    plot_var <- hscn$data %>%
      dplyr::group_by(state, state_AS_phased) %>%
      dplyr:: filter(Min > 0 & Maj > 0) %>%
      dplyr::summarize(mBAF = median(BAF), varBAF = var(BAF)) %>%
      ggplot2::ggplot(ggplot2::aes(x = state, y = varBAF, col = factor(paste0("CN", state), levels = paste0("CN", seq(0, 11, 1))))) +
      ggplot2::geom_text(ggplot2::aes(label = state_AS_phased)) +
      cowplot::theme_cowplot() +
      ggplot2::xlab("State") +
      ggplot2::ylab("BAF variance") +
      ggplot2::scale_color_manual(name = "Copy number \n state",
                                 breaks = paste0("CN", seq(0, 11, 1)),
                                 labels = paste0("CN", seq(0, 11, 1)),
                                 values = scCN_cols(paste0("CN", seq(0, 11, 1))),
                                 drop = FALSE) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(labels = seq(1, max(hscn$data$state), 1),
                                  breaks = seq(1, max(hscn$data$state), 1))
  } else {
    plot_var <- hscn$data %>%
      dplyr:: filter(Min > 0 & Maj > 0) %>%
      dplyr::group_by(state, state_AS_phased) %>%
      dplyr::summarize(mBAF = median(BAF), varBAF = var(BAF)) %>%
      dplyr::group_by(state) %>%
      dplyr::summarise(mBAF = mean(mBAF), ymin = min(varBAF), ymax = max(varBAF), varBAF = mean(varBAF)) %>%
      ggplot2::ggplot(ggplot2::aes(x = state, y = varBAF,
                                   ymin = ymin, ymax = ymax,
                                   col = factor(paste0("CN", state), levels = paste0("CN", seq(0, 11, 1))))) +
      ggplot2::geom_pointrange() +
      cowplot::theme_cowplot() +
      ggplot2::xlab("State") +
      ggplot2::ylab("BAF variance") +
      ggplot2::scale_color_manual(name = "Copy number \n state",
                                  breaks = paste0("CN", seq(0, 11, 1)),
                                  labels = paste0("CN", seq(0, 11, 1)),
                                  values = scCN_cols(paste0("CN", seq(0, 11, 1))),
                                  drop = FALSE) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(labels = seq(1, max(hscn$data$state), 1),
                                  breaks = seq(1, max(hscn$data$state), 1))
  }

  return(plot_var)
}

