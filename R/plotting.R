plottinglist <- function(CNbins){
  #arrange segments in order, generate segment index and reorder CN state factor

  chridx <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))
  CNbins <- CNbins %>%
    dplyr::left_join(., chridx, by = c("chr")) %>%
    dplyr::arrange(cell_id, idx, start) %>%
    dplyr::group_by(cell_id) %>%
    dplyr::mutate(idx = 1:dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>%
    dplyr::mutate(CNs = forcats::fct_reorder(ifelse(is.na(state), NA,
                                                    paste0("CN", state)), state))

  #get breaks - first index of each chromosome
  chrbreaks <- CNbins %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(idx)

  chrlabels <- gtools::mixedsort(unique(CNbins$chr))

  return(list(CNbins = CNbins, chrbreaks = chrbreaks, chrlabels = chrlabels))
}

#' @export
plot_umap <- function(clustering, bycol = NULL){
  ggplot2::ggplot(clustering, ggplot2::aes(x = umap1, y = umap2)) +
    ggrastr::geom_point_rast(ggplot2::aes_string(col = bycol)) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_bw() +
    guides(colour = ggplot2::guide_legend(ncol = 3))
}

#' @export
plot_umap <- function(clustering, bycol = NULL, alphavalue = 0.5){
  ggplot2::ggplot(clustering, ggplot2::aes(x = umap1, y = umap2)) +
    ggrastr::geom_point_rast(ggplot2::aes_string(col = bycol), alpha = alphavalue) +
    ggplot2::xlab("UMAP 1") +
    ggplot2::ylab("UMAP 2") +
    ggplot2::theme_bw() +
    guides(colour = ggplot2::guide_legend(ncol = 3))
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
                         returnlist = FALSE){
  if (is.null(cellid)){
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)){
    message(paste0("Filtering for chromosome: ", chrfilt))
    CNbins <- dplyr::filter(CNbins, chr == chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(.)

  gCN <- pl$CNbins %>%
    dplyr::mutate(state = paste0(state)) %>%
    ggplot2::ggplot(ggplot2::aes(x = idxs, y = copy)) +
    ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
    ggplot2::scale_color_manual(name = "Allele Specific CN",
                                breaks = names(statecolpal),
                                labels = names(statecolpal),
                                values = statecolpal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(breaks = pl$chrbreaks, labels = pl$chrlabels) +
    ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN)) +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Copy Number") +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "x") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha=1, size = 2))) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")

  return(gCN)
}


#' @export
plotCNprofileBAF <- function(CNbins,
                          cellid = NULL,
                          chrfilt = NULL,
                          pointsize = 1,
                          alphaval = 0.9,
                          maxCN = 10,
                          cellidx = 1,
                          BAFcol = "state_phase",
                          statecol = "state",
                          returnlist = FALSE){
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
    CNbins <- dplyr::filter(CNbins, chr == chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(.)

  gBAF <- pl$CNbins %>%
    dplyr::mutate(state_min = paste0(state_min)) %>%
    ggplot2::ggplot(ggplot2::aes(x = idxs, y = BAF)) +
    ggplot2::geom_point(ggplot2::aes_string(col = BAFcol), size = pointsize, alpha = alphaval) +
    ggplot2::scale_color_manual(name = "CN",
                                breaks = names(BAFcolpal),
                                labels = names(BAFcolpal),
                                values = BAFcolpal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(breaks = pl$chrbreaks, labels = pl$chrlabels) +
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
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha=1, size = 2)))

  gCN <- pl$CNbins %>%
    dplyr::mutate(state = paste0(state)) %>%
    dplyr::mutate(state_min = paste0(state_min)) %>%
    ggplot2::ggplot(ggplot2::aes(x = idxs, y = copy)) +
    ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval) +
    ggplot2::scale_color_manual(name = "Allele Specific CN",
                                breaks = names(statecolpal),
                                labels = names(statecolpal),
                                values = statecolpal) +
    ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::scale_x_discrete(breaks = pl$chrbreaks, labels = pl$chrlabels) +
    ggplot2::scale_y_continuous(breaks = seq(0, maxCN, 2), limits = c(0, maxCN)) +
    ggplot2::xlab("Chromosome") +
    ggplot2::ylab("Copy Number") +
    cowplot::theme_cowplot() +
    cowplot::background_grid(major = "x") +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha=1, size = 2))) +
    ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = "bottom")

  g <- cowplot::plot_grid(gBAF, gCN, align = "v", ncol = 1, rel_heights = c(1, 1.2))

  if (returnlist == TRUE){
    p <- list(CN = gCN, BAF = gBAF, both = g)
  } else {
    p <- g
  }

  return(p)
}

#' @export
plotCNBAF <- function(CNBAF, nfilt = 10^5, plottitle = "5Mb", pointsize = 0.1){
  maj <- seq(1, 11, 1)
  min <- seq(0, 11, 1)
  ASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS = paste0(state-min, "|", min)) %>%
    dplyr::mutate(Maj = state-min, Min = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(Maj >= 0, Min >= 0)

  CNBAF <- CNBAF %>%
    dplyr::filter(state < 10, copy < 10)

  if (nfilt < length(CNBAF$chr)){
    CNBAF <- CNBAF %>%
      dplyr::sample_n(nfilt)
  }

  plottitle <- paste0((CNBAF$end[1] - CNBAF$start[1] + 1) / 1e6, "Mb bins")

  g <- CNBAF %>%
    ggplot2::ggplot(ggplot2::aes(y = BAF, x = copy, col = paste0("CN", state))) +
    ggplot2::geom_point(size = pointsize, alpha = 0.2) +
    ggplot2::xlab("Corrected read counts") +
    cowplot::theme_cowplot() +
    ggplot2::ggtitle(plottitle) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size=5))) +
    ggplot2::geom_text(data = ASstates, ggplot2::aes(x = state, y = cBAF, label = state_AS), col = "black") +
    ggplot2::scale_color_manual(name = "Copy number \n state",
                       breaks = paste0("CN", seq(0, max(CNBAF$state, na.rm = TRUE), 1)),
                       labels = seq(0, max(CNBAF$state, na.rm = TRUE), 1),
                       values = scCN_cols(paste0("CN", seq(0, max(CNBAF$state, na.rm = TRUE), 1)))) +
    ggplot2::scale_x_continuous(breaks = seq(0, 10, 1), labels = seq(0, 10, 1), limits = c(0, 10))
  return(g)
}

#' @export
plotBAFperstate <- function(alleleCN, minpts = 250, maxstate = 6){

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
    dplyr::filter(state <= maxstate, n > minpts)
  allASstates <- allASstates %>%
    dplyr::filter(state_AS_phased %in% unique(forplot$state_AS_phased))

  g <- forplot %>%
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_reorder(state_AS_phased, state),
                                 y = BAF,
                                 fill = paste0("CN", state))) +
    ggplot2::geom_violin(scale = "width", col = "white") +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, col = "white") +
    ggplot2::scale_fill_manual(name = "Copy number \n state",
                                breaks = paste0("CN", seq(0, max(alleleCN$state, na.rm = TRUE), 1)),
                                labels = seq(0, max(alleleCN$state, na.rm = TRUE), 1),
                                values = scCN_cols(paste0("CN", seq(0, max(alleleCN$state, na.rm = TRUE), 1)))) +
    cowplot::theme_cowplot() +
    ggplot2::geom_crossbar(data = allASstates, ggplot2::aes(y = cBAF, ymin = cBAF, ymax = cBAF)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(g)
}



