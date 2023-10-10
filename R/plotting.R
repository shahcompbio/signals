removexaxis <- ggplot2::theme(axis.line.x=ggplot2::element_blank(),
                     axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank())

removeyaxis <- ggplot2::theme(axis.line.y=ggplot2::element_blank(),
                     axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank())

plottinglist <- function(CNbins, 
                         xaxis_order = "genome_position", 
                         maxCN = 20, 
                         tickwidth = 50, 
                         chrstart = NULL, 
                         positionticks = FALSE,
                         chrend = NULL) {
  # arrange segments in order, generate segment index and reorder CN state factor

  if (!xaxis_order %in% c("bin", "genome_position")) {
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }
  
  if (nrow(CNbins) == 0){
    stop("Data is empty!")
  }
  
  binsize <- CNbins$end[1] - CNbins$start[1] + 1
  tickwidth <- (tickwidth * 1e6) / binsize
  
  if (xaxis_order == "bin") {

    bins <- getBins(binsize = binsize) %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::mutate(idx = 1:dplyr::n())

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
        paste0("CN", state)
      ), state))

    # get breaks - first index of each chromosome
    chrbreaks <- CNbins %>%
      dplyr::group_by(chr) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(idx)

    # get ticks - median bin of each chromosome
    chrticks <- CNbins %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::group_by(chr) %>%
      dplyr::summarise(idx = round(median(idx))) %>%
      dplyr::pull(idx)

    chrlabels <- gtools::mixedsort(unique(CNbins$chr))

    minidx <- min(CNbins$idx)
    maxidx <- max(CNbins$idx)
  } else {

    bins <- getBins(binsize = binsize) %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::mutate(idx = 1:dplyr::n())
    
    if(dim(bins)[1] < tickwidth){
      #ensure width of ticks is larget than number of bins
      tickwidth <- tickwidth / 2
    }

    CNbins <- dplyr::full_join(bins, CNbins, by = c("chr", "start", "end")) %>%
      dplyr::filter(!is.na(copy)) %>%
      dplyr::filter(!is.na(state)) %>%
      dplyr::mutate(copy = ifelse(copy > maxCN, maxCN, copy)) %>%
      #dplyr::mutate(state = ifelse(state > maxCN, maxCN, state)) %>%
      dplyr::mutate(idxs = forcats::fct_reorder(factor(idx), idx)) %>%
      dplyr::mutate(CNs = forcats::fct_reorder(ifelse(is.na(state), NA,
        paste0("CN", state)
      ), state))

    # get breaks - first index of each chromosome
    chrbreaks <- bins %>%
      dplyr::filter(chr %in% unique(CNbins$chr)) %>%
      dplyr::group_by(chr) %>%
      dplyr::filter(dplyr::row_number() == 1) %>%
      dplyr::pull(idx)

    # get ticks - median bin of each chromosome
    if (length(unique(CNbins$chr)) == 1){
      chrticks <- seq(tickwidth, dim(bins)[1], tickwidth)
      chrlabels <- paste0((chrticks * binsize) / 1e6)
    } else {
      chrticks <- bins %>%
        dplyr::filter(chr %in% unique(CNbins$chr)) %>%
        dplyr::group_by(chr) %>%
        dplyr::summarise(idx = round(median(idx))) %>%
        dplyr::pull(idx)
      chrlabels <- gtools::mixedsort(unique(CNbins$chr))
    }
    
    if (length(unique(CNbins$chr)) == 1 & (!is.null(chrstart) | !is.null(chrstart))){
      idxstart <- ifelse(is.null(chrstart), min(CNbins$idx), (chrstart * 1e6) / binsize)
      idxend <- ifelse(is.null(chrend), max(CNbins$idx), (chrend * 1e6) / binsize)
      CNbins <- CNbins %>% dplyr::filter(idx >= idxstart)
      CNbins <- CNbins %>% dplyr::filter(idx <= idxend)
      chrticks <- seq(idxstart, idxend, tickwidth)
      chrlabels <- paste0((chrticks * binsize) / 1e6)
      minidx <- ifelse(is.null(chrstart), min(bins$idx), idxstart)
      maxidx <- ifelse(is.null(chrend), max(bins$idx), idxend)
    } else{
      minidx <- min(bins$idx)
      maxidx <- max(bins$idx)
    }
    
    getticks <- function(mychr){
      CNbins_filt <- dplyr::filter(CNbins, chr == mychr)
      bins_filt <- dplyr::filter(bins, chr == mychr)
      idxstart <- min(CNbins_filt$idx)
      idxend <- max(CNbins_filt$idx)
      chrticks <- seq(tickwidth + bins_filt$idx[1], tail(bins_filt$idx,1), tickwidth)
      chrlabels <- paste0(((chrticks - bins_filt$idx[1]) * binsize) / 1e6)
      return(list(chrlabels = chrlabels, chrticks = chrticks))
    }
    
    if (positionticks){
      idxstart <- ifelse(is.null(chrstart), min(CNbins$idx), (chrstart * 1e6) / binsize)
      idxend <- ifelse(is.null(chrend), max(CNbins$idx), (chrend * 1e6) / binsize)
      CNbins <- CNbins %>% dplyr::filter(idx >= idxstart)
      CNbins <- CNbins %>% dplyr::filter(idx <= idxend)
      dat <- lapply(unique(CNbins$chr), getticks)
      chrticks <- unlist(lapply(dat, function(x) x["chrticks"][[1]]))
      chrlabels <- unlist(lapply(dat, function(x) x["chrlabels"][[1]]))
      minidx <- ifelse(is.null(chrstart), min(bins$idx), idxstart)
      maxidx <- ifelse(is.null(chrend), max(bins$idx), idxend)
    }

  }

  return(list(CNbins = CNbins, bins = bins, chrbreaks = chrbreaks, chrticks = chrticks, chrlabels = chrlabels, minidx = minidx, maxidx = maxidx))
}

plottinglistSV <- function(breakpoints, binsize = 0.5e6, chrfilt = NULL, chrstart = NULL, chrend = NULL) {
  breakpoints$chromosome_1 <- as.character(breakpoints$chromosome_1)
  breakpoints$chromosome_2 <- as.character(breakpoints$chromosome_2)

  if (is.null(chrfilt)) {
    bins <- getBins(binsize = binsize) %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(idx = 1:dplyr::n()) %>%
      dplyr::select(chr, start, idx) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(maxidx = max(idx), minidx = min(idx)) %>%
      dplyr::ungroup()
  } else {
    bins <- getBins(binsize = binsize) %>%
      dplyr::filter(chr %in% chrfilt) %>%
      dplyr::filter(chr != "Y") %>%
      dplyr::mutate(idx = 1:dplyr::n()) %>%
      dplyr::select(chr, start, idx) %>%
      dplyr::group_by(chr) %>%
      dplyr::mutate(maxidx = max(idx), minidx = min(idx)) %>%
      dplyr::ungroup()
    breakpoints <- breakpoints %>%
      dplyr::filter((chromosome_1 %in% chrfilt) | (chromosome_2 %in% chrfilt))
  }

  breakpoints <- breakpoints %>%
    dplyr::mutate(
      position_1 = ifelse(position_1 < position_2, 
                          binsize * floor(position_1 / binsize) + 1,
                          binsize * ceiling(position_1 / binsize) + 1),
      position_2 = ifelse(position_2 < position_1, 
                          binsize * floor(position_2 / binsize) + 1,
                          binsize * ceiling(position_2 / binsize) + 1)
    ) %>%
    dplyr::mutate(position_2 = ifelse(abs(position_1- position_2) <= binsize, position_1, position_2)) %>% 
    dplyr::left_join(bins %>% dplyr::rename(chromosome_1 = chr, position_1 = start, idx_1 = idx, maxidx_1 = maxidx, minidx_1 = minidx), by = c("chromosome_1", "position_1")) %>%
    dplyr::left_join(bins %>% dplyr::rename(chromosome_2 = chr, position_2 = start, idx_2 = idx, maxidx_2 = maxidx, minidx_2 = minidx), by = c("chromosome_2", "position_2"))

  breakpoints <- breakpoints %>%
    dplyr::distinct(chromosome_1, position_1, chromosome_2, position_2, type, rearrangement_type, idx_1, idx_2, maxidx_1, minidx_1, maxidx_2, minidx_2)

  # get breaks - first index of each chromosome
  chrbreaks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::pull(idx)

  # get ticks - median bin of each chromosome
  chrticks <- bins %>%
    dplyr::filter(chr %in% unique(CNbins$chr)) %>%
    dplyr::group_by(chr) %>%
    dplyr::summarise(idx = round(median(idx))) %>%
    dplyr::pull(idx)

  chrlabels <- gtools::mixedsort(unique(bins$chr))
  if (length(chrfilt) == 1 & (!is.null(chrstart) | !is.null(chrstart))){
    minidx <- ifelse(is.null(chrstart), min(bins$idx), (chrstart * 1e6) / binsize)
    maxidx <- ifelse(is.null(chrend), max(bins$idx), (chrend * 1e6) / binsize)
  } else{
    minidx <- min(bins$idx)
    maxidx <- max(bins$idx)
  }
  return(list(breakpoints = breakpoints, bins = bins, chrbreaks = chrbreaks, chrticks = chrticks, chrlabels = chrlabels, minidx = minidx, maxidx = maxidx))
}

CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1, 1)), substring(c, 2),
    sep = "", collapse = " "
  )
}

#' @export
plotSV <- function(breakpoints,
                   chrfilt = NULL,
                   curvature = -0.5,
                   returnlist = FALSE,
                   ylims = c(0, 2),
                   legend.position = "bottom",
                   ...) {
  pl <- plottinglistSV(breakpoints, chrfilt = chrfilt)

  pl$breakpoints <- pl$breakpoints %>%
    dplyr::mutate(curve = ifelse(abs(idx_1 - idx_2) < 5, FALSE, TRUE))

  pl$breakpoints$rearrangement_type <- unlist(lapply(pl$breakpoints$rearrangement_type, CapStr))

  curve_data <- pl$breakpoints %>% dplyr::filter(curve == TRUE)
  line_data <- pl$breakpoints %>% dplyr::filter(curve == FALSE)

  gSV <- pl$bins %>%
    ggplot2::ggplot(ggplot2::aes(x = idx, y = 1)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) +
    xlab("Chromosome") +
    cowplot::theme_cowplot(...) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = legend.position
    ) +
    ylab("SV") +
    ggplot2::ylim(ylims)

  if (dim(curve_data)[1] > 0) {
    gSV <- gSV + ggplot2::geom_curve(data = curve_data, ggplot2::aes(x = idx_1, xend = idx_2, y = 1, yend = 1.0001, col = rearrangement_type), curvature = curvature) +
      ggplot2::labs(col = "Rearrangement") +
      ggplot2::scale_color_manual(
        breaks = names(SV_colors),
        values = as.vector(SV_colors)
      )
  }

  if (dim(line_data)[1] > 0) {
    line_data <- line_data %>%
      dplyr::mutate(y = ifelse(rearrangement_type == "Foldback", 1, 1)) %>%
      dplyr::mutate(yend = ifelse(rearrangement_type == "Foldback", min(ylims), max(ylims)))

    gSV <- gSV + ggplot2::geom_segment(data = line_data, ggplot2::aes(x = idx_1, xend = idx_1 + 0.001, y = y, yend = yend, col = rearrangement_type)) +
      ggplot2::labs(col = "Rearrangement") +
      ggplot2::scale_color_manual(
        breaks = names(SV_colors),
        values = as.vector(SV_colors)
      )
  }

  if (returnlist == TRUE) {
    p <- list(SV = gSV, plist = pl)
  } else {
    p <- gSV
  }

  return(p)
}


#' @export
plotSV2 <- function(breakpoints,
                    chrfilt = NULL,
                    curvature = -0.5,
                    returnlist = FALSE,
                    ylims = c(0, 2),
                    legend.position = "bottom",
                    svwidth = 1.0,
                    chrstart = NULL,
                    chrend = NULL,
                    ...) {
  pl <- plottinglistSV(breakpoints, chrfilt = chrfilt, chrstart = chrstart, chrend = chrend)

  pl$breakpoints <- pl$breakpoints %>%
    dplyr::mutate(curve = ifelse((abs(idx_1 - idx_2) < 2) & (chromosome_1 == chromosome_2), FALSE, TRUE))

  pl$breakpoints$rearrangement_type <- unlist(lapply(pl$breakpoints$rearrangement_type, CapStr))

  if (dim(pl$breakpoints)[1] > 0) {
    bezdf <- get_bezier_df_2(pl$breakpoints %>% dplyr::filter(curve == TRUE))
  }

  gSV <- pl$bins %>%
    ggplot2::ggplot(ggplot2::aes(x = idx, y = 1)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
    ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx)) +
    xlab("Chromosome") +
    cowplot::theme_cowplot(...) +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = legend.position
    ) +
    ylab("SV") +
    ggplot2::ylim(ylims)

  if (dim(pl$breakpoints)[1] > 0) {
    gSV <- gSV +
      ggforce::geom_bezier(ggplot2::aes(x = idx, y = yidx, group = id, col = rearrangement_type),
        alpha = 0.5,
        size = svwidth,
        data = bezdf
      ) +
      ggplot2::labs(col = "Rearrangement") +
      ggplot2::scale_color_manual(
        breaks = names(SV_colors),
        values = as.vector(SV_colors)
      )
  }

  line_data <- pl$breakpoints %>% dplyr::filter(curve == FALSE)

  if (dim(line_data)[1] > 0) {
    line_data <- line_data %>%
      dplyr::mutate(y = ifelse(rearrangement_type == "Foldback", 1, 1)) %>%
      dplyr::mutate(yend = ifelse(rearrangement_type == "Foldback", min(ylims), max(ylims)))

    gSV <- gSV + ggplot2::geom_segment(data = line_data, ggplot2::aes(x = idx_1, xend = idx_1 + 0.001, y = y, yend = yend, col = rearrangement_type), size = svwidth) +
      ggplot2::labs(col = "Rearrangement") +
      ggplot2::scale_color_manual(
        breaks = names(SV_colors),
        values = as.vector(SV_colors)
      )
  }

  if (returnlist == TRUE) {
    p <- list(SV = gSV, plist = pl)
  } else {
    p <- gSV
  }

  return(p)
}

get_gene_idx <- function(mygenes, chr = NULL, binsize = 0.5e6) {
  gene_df <- gene_locations %>%
    dplyr::filter(ensembl_gene_symbol %in% mygenes) %>%
    dplyr::mutate(
      start = binsize * floor(start / binsize) + 1,
      end = binsize * ceiling(start / binsize)
    )
  bins <- getBins(binsize = binsize, chromosomes = chr) %>% dplyr::mutate(idx = 1:dplyr::n())
  gene_bin <- dplyr::left_join(gene_df, bins)
  return(gene_bin)
}

#' @export
plot_umap <- function(clustering, bycol = NULL, alphavalue = 0.5, raster = FALSE) {
  if (raster == TRUE) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
        call. = FALSE
      )
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
    function(x) tanh(0.075 * x),
    function(x) {
      nnan <- 1
      while (nnan > 0){
        suppressWarnings(y <- atanh(x) / 0.075)
        nnan <- sum(is.nan(y))
        x[is.nan(y)] <- x[is.nan(y)] - 0.000001
      }
      return(y)
      }
    
  )
}

#create_squashy_trans

get_bezier_df_2 <- function(sv) {
  set.seed(123)

  suppressWarnings({
    sv_ <- sv %>%
      dplyr::filter(chromosome_1 != "X") %>%
      dplyr::mutate(outer = is.na(idx_2)) %>%
      dplyr::mutate(idx_2 = ifelse(outer & (as.numeric(chromosome_1) < as.numeric(chromosome_2)), maxidx_1, idx_2)) %>%
      dplyr::mutate(idx_2 = ifelse(outer & (as.numeric(chromosome_1) > as.numeric(chromosome_2)), minidx_1, idx_2)) %>%
      dplyr::mutate(idx_2 = ifelse(is.na(idx_2), maxidx_1, idx_2)) %>%
      dplyr::mutate(idx_2 = ifelse(is.na(idx_1), maxidx_1, idx_2)) %>%
      dplyr::mutate(idx_3 = idx_2) %>%
      dplyr::mutate(idx_2 = (idx_1 + idx_3) / 2) %>%
      dplyr::mutate(yidx_1 = 1, yidx_2 = 2, yidx_3 = ifelse(outer, 2, 1)) %>%
      dplyr::mutate(yidx_2 = ifelse(rearrangement_type == "Foldback", 0, 2)) %>%
      dplyr::select(-maxidx_1, -maxidx_2, -minidx_1, -minidx_2)
  })

  x1 <- sv_ %>%
    dplyr::select(rearrangement_type, idx_1, idx_2, idx_3, chromosome_1, chromosome_2, position_1, position_2) %>%
    tidyr::pivot_longer(dplyr::starts_with("idx"), names_to = "name1", values_to = "idx") %>%
    dplyr::arrange(chromosome_1, chromosome_2, position_1, position_2, name1) %>%
    tidyr::separate(name1, c("i", "bezidx"), sep = "_") %>%
    dplyr::select(-i)

  x2 <- sv_ %>%
    dplyr::select(rearrangement_type, yidx_1, yidx_2, yidx_3, chromosome_1, position_1, position_2) %>%
    tidyr::pivot_longer(dplyr::starts_with("yidx"), names_to = "name2", values_to = "yidx") %>%
    dplyr::arrange(chromosome_1, position_1, name2) %>%
    tidyr::separate(name2, c("i", "bezidx"), sep = "_") %>%
    dplyr::select(-i)

  bez <- dplyr::left_join(x1, x2,
    by = c(
      "rearrangement_type", "chromosome_1",
      "position_1", "position_2", "bezidx"
    )
  ) %>%
    dplyr::mutate(id = paste(chromosome_1, position_1, position_2,
      rearrangement_type,
      sep = "_"
    )) %>%
    dplyr::distinct(.)

  return(bez)
}

get_bezier_df <- function(sv, cn, maxCN, homolog = FALSE) {
  set.seed(123)

  if (homolog == TRUE) {
    cn$CNbins <- cn$CNbins %>%
      dplyr::mutate(copy = pmax(Acopy, Bcopy))
  }

  maxidx <- max(cn$CNbins$idx, na.rm = TRUE)
  minidx <- min(cn$CNbins$idx, na.rm = TRUE)
  idxrange <- 0.05 * maxidx
  
  cn$CNbins <- cn$CNbins %>% 
    tidyr::fill( c("state", "copy"), .direction = "downup")

  svcn1 <- dplyr::left_join(sv$breakpoints, cn$CNbins %>%
    dplyr::rename(chromosome_1 = chr, position_1 = start, copy_1 = copy) %>%
    dplyr::select(chromosome_1, position_1, copy_1))
  svcn <- dplyr::inner_join(svcn1, cn$CNbins %>%
    dplyr::rename(chromosome_2 = chr, position_2 = start, copy_2 = copy) %>%
    dplyr::select(chromosome_2, position_2, copy_2)) %>%
    dplyr::mutate(copy_1 = ifelse(is.na(copy_1), 2, copy_1)) %>%
    dplyr::mutate(copy_2 = ifelse(is.na(copy_2), 2, copy_2)) %>%
    dplyr::distinct(.) %>%
    na.omit(.) %>%
    dplyr::rename(idx_3 = idx_2, copy_3 = copy_2) %>%
    mutate(copy_1 = ifelse(position_1 == position_2, 0, copy_1))

  x1 <- svcn %>%
    dplyr::select(rearrangement_type, idx_1, idx_3, chromosome_1, chromosome_2, position_1, position_2) %>%
    dplyr::mutate(minidx = pmin(idx_1, idx_3)) %>%
    dplyr::mutate(sr = ifelse(rearrangement_type == "foldback", 0,
      sample(c(-1, 1), dplyr::n(), replace = TRUE) * idxrange
    )) %>%
    dplyr::mutate(sr = dplyr::case_when(
      rearrangement_type == "foldback" ~ 0,
      idx_3 > idx_1 ~ 1,
      idx_3 < idx_1 ~ -1,
      idx_3 == idx_1 ~ 0,
    )) %>%
    dplyr::mutate(idx_2 = minidx + sr) %>%
    tidyr::pivot_longer(dplyr::starts_with("idx"), names_to = "name1", values_to = "idx") %>%
    dplyr::arrange(chromosome_1, position_1, position_2, name1) %>%
    tidyr::separate(name1, c("i", "bezidx"), sep = "_") %>%
    dplyr::select(-sr, -i)

  x2 <- svcn %>%
    dplyr::select(rearrangement_type, copy_1, copy_3, chromosome_1, position_1, position_2) %>%
    dplyr::mutate(mincopy = pmin(copy_1, copy_3)) %>%
    dplyr::mutate(sr1 = sample(c(0.25, 0.75), dplyr::n(), replace = TRUE)) %>%
    dplyr::mutate(sr2 = sample(c(-1, 1), dplyr::n(), replace = TRUE)) %>%
    dplyr::mutate(copydiff = abs(copy_1 - copy_3)) %>%
    dplyr::mutate(copy_2 = ifelse(copydiff > 3, mincopy + sr1 * copydiff, mincopy + sr2 * 2.5)) %>%
    dplyr::mutate(copy_2 = ifelse(copy_2 > maxCN, maxCN, copy_2)) %>%
    dplyr::mutate(copy_2 = ifelse(copy_2 < 0, 1, copy_2)) %>%
    dplyr::mutate(copy_2 = ifelse(rearrangement_type == "foldback", copy_3, copy_2)) %>%
    dplyr::select(-copydiff) %>%
    tidyr::pivot_longer(dplyr::starts_with("copy"), names_to = "name2", values_to = "copy") %>%
    dplyr::arrange(chromosome_1, position_1, name2) %>%
    tidyr::separate(name2, c("i", "bezidx"), sep = "_") %>%
    dplyr::select(-sr1, -sr2, -i)

  bez <- dplyr::left_join(x1, x2,
    by = c(
      "rearrangement_type", "chromosome_1",
      "position_1", "position_2", "bezidx"
    )
  ) %>%
    dplyr::mutate(id = paste(chromosome_1, position_1, position_2,
      rearrangement_type,
      sep = "_"
    )) %>%
    dplyr::distinct(.) %>%
    dplyr::mutate(idx = ifelse(idx > maxidx, maxidx, idx)) %>%
    dplyr::mutate(idx = ifelse(idx < minidx, minidx, idx))

  return(bez)
}

#' Plot a single cell copy number profile
#'
#' @param CNbins Single cell copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`
#' @param cellid Which cell to plot, if no cell is specific will plot the first cell in the dataframe
#' @param chrfilt Vector of chromosomes to plot, if NULL (default) will plot all chromosomes
#' @param pointsize The point size in the plot
#' @param alphaval Alpha value of points
#' @param maxCN The maximum on the y axis, if any points are above this value they will be winsorized rather than removed
#' @param cellidx idx of cell to plot if cellid = NULL
#' @param statecol The colour mapping, default is to map colours to the `state` column
#' @param returnlist Return a list rather than the ggplot object
#' @param raster use ggrastr or not, default = FALSE
#' @param y_axis_trans What transformation to use on the y-axis, default is identity, the other option is "squashy" which uses a tanh transformation
#' @param xaxis_order Default is "genome_position"
#' @param legend.position Where to place the legend, default is "bottom"
#' @param annotateregions Dataframe with chr start and end positions to annotate, will draw a dashed vertical line at this position
#' @param annotateregions_linetype linetype for region annotation, default = 2 (dashed)
#' @param SV Default is NULL. If a dataframe with structural variant position is passed it will add rearrangement links between bins.
#' @param SVcol Default is TRUE. Colour SVs or not
#' @param svalpha the alpha scaling of the SV lines, default = 0.5
#' @param genes vector of genes to annotate, will add a dashed vertical line and label
#' @param tickwidth Spacing of ticks (in Mb) when only 1 chromosome is plotted
#' @param svwidth Width of SV width curves, default = 1.0
#' @param adj adjustment for gene labels
#' @param chrstart Start of region (in Mb) when plotting a single chromosome
#' @param chrend End of region (in Mb) when plotting a single chromosome
#' @param shape shape for plotting, default = 16
#' @param positionticks set to TRUE to use position ticks rather than chromosome ticks
#' @param genome genome to use, default = "hg19" (only used for ideogram)
#' @param ideogram plot ideogram at the top, default = TRUE
#' @param ideogram_height height of the ideogram
#' 
#' @return ggplot2 plot
#'
#' @examples
#'
#' plotCNprofile(CNbins)
#' @md
#' @export
plotCNprofile <- function(CNbins,
                          cellid = NULL,
                          chrfilt = NULL,
                          pointsize = 1,
                          alphaval = 0.6,
                          maxCN = 10,
                          cellidx = 1,
                          statecol = "state",
                          returnlist = FALSE,
                          raster = FALSE,
                          genome = "hg19",
                          y_axis_trans = "identity",
                          xaxis_order = "genome_position",
                          legend.position = "bottom",
                          annotateregions = NULL,
                          annotateregions_linetype = 2,
                          SV = NULL,
                          SVcol = TRUE,
                          svalpha = 0.5,
                          svwidth = 1.0,
                          adj = 0.03,
                          genes = NULL, 
                          tickwidth = 50,
                          chrstart = NULL,
                          chrend = NULL,
                          shape = 16,
                          positionticks = FALSE,
                          ideogram = FALSE,
                          overwrite_color = NULL,
                          ...) {
  if (!xaxis_order %in% c("bin", "genome_position")) {
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }
  
  if (is.null(cellid)) {
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  if (y_axis_trans == "squashy") {
    ybreaks <- c(0, 2, 5, 10, maxCN)
  } else {
    ybreaks <- seq(0, maxCN, 2)
  }
  
  if (length(chrfilt) == 1){
    xlab <- paste0('Chr. ', chrfilt, " (Mb)")
  } else{
    xlab <- 'Chromosome'
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)) {
    message(paste0("Filtering for chromosomes: ", paste0(chrfilt, collapse = ",")))
    CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN, positionticks = positionticks,
                 tickwidth = tickwidth, chrstart = chrstart, chrend = chrend)
  
  if (ideogram == TRUE){
    miny <- -0.5
  } else{
    miny <- 0
  }
  
  if (raster == TRUE) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }
    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = shape) +
      ggplot2::scale_color_manual(
        name = "Copy number",
        breaks = names(statecolpal),
        labels = names(statecolpal),
        values = statecolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(miny, maxCN), trans = y_axis_trans) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot(...) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = 6, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3, shape = 15)
      )) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
  } else {
    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = 16) +
      ggplot2::scale_color_manual(
        name = "Allele Specific CN",
        breaks = names(statecolpal),
        labels = names(statecolpal),
        values = statecolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(miny, maxCN), trans = y_axis_trans) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot(...) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = 6, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3, shape = 15)
      )) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
  }

  if (!is.null(genes)) {
    yplace <- maxCN
    if (ideogram == TRUE){
      yplace <- yplace - 2
    }
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    gene_idx <- get_gene_idx(genes, chr = chrfilt, binsize = binsize)
    npoints <- dim(pl$CNbins)[1]
    gCN <- gCN +
      ggplot2::geom_vline(data = gene_idx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3) +
      ggrepel::geom_text_repel(data = gene_idx, ggplot2::aes(x = idx - npoints * adj, y = maxCN, label = ensembl_gene_symbol), col = "black", alpha = 0.75)
  }

  if (!is.null(annotateregions)) {
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    annotateregions <- dplyr::mutate(annotateregions, start = round(start / binsize) * binsize + 1)
    datidx <- dplyr::inner_join(annotateregions, pl$bins %>% dplyr::select(chr, start, idx)) %>% dplyr::distinct(.)
    datidx <- dplyr::mutate(datidx, idx = ifelse(idx > pl$maxidx, pl$maxidx, idx))
    gCN <- gCN +
      ggplot2::geom_vline(data = datidx, ggplot2::aes(xintercept = idx), lty = annotateregions_linetype, size = 0.3, alpha = 0.5)
  }
  
  if (!is.null(SV) & !is.null(chrfilt)){
    SV <- dplyr::filter(SV, chromosome_1 %in% chrfilt & chromosome_2 %in% chrfilt)
  }

  if (!is.null(SV) && nrow(SV) > 0) {
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    svpl <- plottinglistSV(SV, chrfilt = chrfilt, binsize = binsize)
    pl$CNbins <- dplyr::left_join(getBins(binsize = binsize), pl$CNbins)
    bezdf <- get_bezier_df(svpl, pl, maxCN)
    bezdf <- bezdf %>%
      dplyr::mutate(samebin = (position_1 == position_2) | 
                      rearrangement_type == "foldback")
    bezdf$rearrangement_type = unlist(lapply(bezdf$rearrangement_type, CapStr))
    
    rearrangement_types <- unique(bezdf %>% pull(rearrangement_type))
    
    if (SVcol == TRUE){ #different colours for each SV
      for (r in rearrangement_types){
        gCN <- gCN +
          ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                               alpha = 0.8,
                               size = svwidth,
                               col = as.vector(SV_colors[r]),
                               data = bezdf %>% dplyr::filter(rearrangement_type == r & samebin == TRUE)
          ) +
          ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                               alpha = svalpha,
                               size = svwidth,
                               col = as.vector(SV_colors[r]),
                               data = bezdf %>% dplyr::filter(rearrangement_type == r & samebin == FALSE)
          )
      }
    } else { # grey colour for arcs, brown color for lines (SVs with start end in the same bins)
      gCN <- gCN +
        ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                             alpha = 0.8,
                             size = svwidth,
                             col = as.vector(SV_colors["Inversion"]),
                             data = bezdf %>% dplyr::filter(samebin == TRUE)
        ) +
        ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
                             alpha = svalpha,
                             size = svwidth,
                             col = "grey30",
                             data = bezdf %>% dplyr::filter(samebin == FALSE)
        )
    }
  }
  
  if (!is.null(overwrite_color)){
    gCN <- gCN + 
      ggplot2::scale_color_manual(values = overwrite_color) +
      ggplot2::theme(legend.position = "none")
  }
  
  if (ideogram == TRUE){
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    ideogram_dat <- cytoband_map[[genome]]
    names(ideogram_dat) <- c("chr", "start", "end", "band", "colval")
    ideogram_dat <- ideogram_dat %>% 
      dplyr::mutate(chr = stringr::str_remove(chr, "chr")) %>% 
      dplyr::mutate(start = round(start / binsize) * binsize + 1, 
                    end = round(end / binsize) * binsize + 1)
    
    #create a dataframe that has the index of the start and end position
    cnbin_idx_start <- pl$bins %>% 
      dplyr::select(chr, start, idx) %>% 
      dplyr::rename(idx_start = idx)
    cnbin_idx_end <-  pl$bins %>% 
      dplyr::select(chr, start, idx) %>% 
      dplyr::rename(end = start) %>% 
      dplyr::rename(idx_end = idx)
    ideogram_dat <- dplyr::inner_join(ideogram_dat,
                                       cnbin_idx_start, by = c("chr", "start")) %>% 
      dplyr::inner_join(cnbin_idx_end, by = c("chr", "end"))
    
    gCN <- gCN +
      ggplot2::geom_rect(data = ideogram_dat,
                         ggplot2::aes(xmin = idx_start, 
                                      y = NULL,
                                      x = NULL,
                             xmax = idx_end, 
                             ymin = -0.5, 
                             ymax = -0.15, fill = colval)) +
      ggplot2::scale_fill_manual(values = cyto_colors) +
      ggplot2::theme(legend.position = "none")
      
  }

  if (returnlist == TRUE) {
    gCN <- list(CN = gCN, plist = pl)
  }

  return(gCN)
}

plotCNprofileBAFhomolog <- function(cn,
                                    cellid = NULL,
                                    chrfilt = NULL,
                                    pointsize = 1,
                                    alphaval = 0.75,
                                    maxCN = 10,
                                    cellidx = 1,
                                    returnlist = FALSE,
                                    raster = FALSE,
                                    y_axis_trans = "identity",
                                    xaxis_order = "genome_position",
                                    legend.position = "bottom",
                                    genes = NULL,
                                    annotateregions = NULL,
                                    annotateregions_linetype = 2,
                                    homolog = FALSE,
                                    SV = NULL,
                                    adj = 0.03,
                                    svalpha = 0.5,
                                    svwidth = 1.0,
                                    shuffle = TRUE,
                                    plotdata = TRUE,
                                    offset = NULL,
                                    linewidth = 0.6,
                                    tickwidth = 50,
                                    chrstart = NULL,
                                    chrend = NULL,
                                    shape = 16,
                                    positionticks = FALSE,
                                    ...) {
  if (!xaxis_order %in% c("bin", "genome_position")) {
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }

  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }
  
  if (length(chrfilt) == 1){
    xlab <- paste0('Chr. ', chrfilt, " (Mb)")
  } else{
    xlab <- 'Chromosome'
  }

  CNbins <- CNbins %>%
    dplyr::mutate(Bcopy = BAF * copy, Acopy = (1 - BAF) * copy)

  if (y_axis_trans == "squashy") {
    ybreaks <- c(0, 2, 5, 10, maxCN)
  } else {
    ybreaks <- seq(0, maxCN, 2)
  }

  if (is.null(cellid)) {
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  if (!"BAF" %in% names(CNbins)) {
    stop("No BAF column in dataframe, first calculate the BAF per bin using combineBAFCN and then callAlleleSpecificCN")
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)) {
    message(paste0("Filtering for chromosomes: ", paste0(chrfilt, collapse = ",")))
    CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    dplyr::mutate(Acopy = ifelse(Acopy > maxCN, maxCN - 0.001, Acopy)) %>%
    dplyr::mutate(Bcopy = ifelse(Bcopy > maxCN, maxCN - 0.001, Bcopy)) %>%
    plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN, positionticks = positionticks,
                 tickwidth = tickwidth, chrstart = chrstart, chrend = chrend)
  pl_for_sv <- pl 

  if (plotdata == TRUE){
    pl$CNbins <- pl$CNbins %>% 
      tidyr::pivot_longer(cols = c("Acopy", "Bcopy"))
    if (shuffle){
      pl$CNbins <- dplyr::sample_frac(pl$CNbins, 1L)
    }
      
    if (raster == TRUE) {
      if (!requireNamespace("ggrastr", quietly = TRUE)) {
        stop("Package \"ggrastr\" needed for this function to work. Please install it.",
          call. = FALSE
        )
      }
  
      gCN <- pl$CNbins %>%
        dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
        dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
        dplyr::mutate(state_min = paste0(state_min)) %>%
        ggplot2::ggplot(ggplot2::aes(x = idx)) +
        ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
        ggrastr::geom_point_rast(ggplot2::aes(y = value, col = name), size = pointsize, alpha = alphaval, shape = shape) +
        ggplot2::scale_color_manual(
          name = "",
          labels = c("Homolog A", "Homolog B"),
          values = as.vector(scCNphase_colors[c("A-Hom", "B-Hom")]),
          drop = FALSE
        ) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          legend.position = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Copy Number") +
        cowplot::theme_cowplot(...) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          ncol = 6, byrow = TRUE,
          override.aes = list(alpha = 1, size = 3, shape = 15)
        )) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    } else {
      gCN <- pl$CNbins %>%
        dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
        dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
        dplyr::mutate(state_min = paste0(state_min)) %>%
        ggplot2::ggplot(ggplot2::aes(x = idx)) +
        ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
        ggplot2::geom_point(ggplot2::aes(y = value, col = name), size = pointsize, alpha = alphaval, shape = shape) +
        ggplot2::scale_color_manual(
          name = "",
          labels = c("Homolog A", "Homolog B"),
          values = as.vector(scCNphase_colors[c("A-Hom", "B-Hom")]),
          drop = FALSE
        ) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          legend.position = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Copy Number") +
        cowplot::theme_cowplot(...) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          ncol = 6, byrow = TRUE,
          override.aes = list(alpha = 1, size = 3, shape = 15)
        )) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    }
  } else{
    
    pl$CNbins <- pl$CNbins %>% 
      dplyr::mutate(rl = data.table::rleid(state_AS_phased)) %>% 
      dplyr::group_by(chr, rl, B, A) %>% 
      dplyr::summarize(idx_start = dplyr::first(idx), idx_end = dplyr::last(idx))
    
    pl$CNbins <- pl$CNbins %>% 
      tidyr::pivot_longer(cols = c("A", "B"))
    
    if (is.null(offset)) {
      offset <- maxCN * 1.1 * 0.03
    }
    
    if (y_axis_trans == "squashy"){
      pl$CNbins <- pl$CNbins %>%
        dplyr::mutate(value = ifelse(name == "A", value + squashy_trans()$transform(2 * value * offset), 
                                     value - squashy_trans()$transform(2 * value * offset)))
    } else{
      pl$CNbins <- pl$CNbins %>%
        dplyr::mutate(value = ifelse(name == "A", value + offset, value - offset))
    }
    
    if (shuffle){
      pl$CNbins <- dplyr::sample_frac(pl$CNbins, 1L)
    }
    
    if (raster == TRUE) {
      if (!requireNamespace("ggrastr", quietly = TRUE)) {
        stop("Package \"ggrastr\" needed for this function to work. Please install it.",
             call. = FALSE
        )
      }
      
      gCN <- pl$CNbins %>%
        ggplot2::ggplot() +
        ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
        ggplot2::geom_linerange(aes(xmin = idx_start, xmax = idx_end, y = value, colour = name), size = linewidth) +
        ggplot2::scale_color_manual(
          name = "",
          labels = c("Homolog A", "Homolog B"),
          values = as.vector(scCNphase_colors[c("A-Hom", "B-Hom")]),
          drop = FALSE
        ) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          legend.position = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0 - offset, maxCN + offset), trans = y_axis_trans) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Copy Number") +
        cowplot::theme_cowplot(...) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          ncol = 6, byrow = TRUE,
          override.aes = list(alpha = 1, size = 3, shape = 15)
        )) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    } else {
      gCN <- pl$CNbins %>%
        ggplot2::ggplot() +
        ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
        ggplot2::geom_linerange(aes(xmin = idx_start, xmax = idx_end, y = value, colour = name), size = linewidth) +
        ggplot2::scale_color_manual(
          name = "",
          labels = c("Homolog A", "Homolog B"),
          values = as.vector(scCNphase_colors[c("A-Hom", "B-Hom")]),
          drop = FALSE
        ) +
        ggplot2::theme(
          axis.title.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          legend.position = "none"
        ) +
        ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
        ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0 - offset, maxCN + offset), trans = y_axis_trans) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab("Copy Number") +
        cowplot::theme_cowplot(...) +
        ggplot2::guides(colour = ggplot2::guide_legend(
          ncol = 6, byrow = TRUE,
          override.aes = list(alpha = 1, size = 3, shape = 15)
        )) +
        ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
    }
    
  }

  if (!is.null(SV)) {
    svpl <- plottinglistSV(SV, chrfilt = chrfilt)
    bezdf <- get_bezier_df(svpl, pl_for_sv, maxCN, homolog = TRUE)
    bezdf <- bezdf %>%
      dplyr::filter((position_1 != position_2) | rearrangement_type == "foldback")
    gCN <- gCN +
      ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
        alpha = 0.5,
        size = svwidth,
        col = as.vector(SV_colors["Foldback"]),
        data = bezdf %>% dplyr::filter(rearrangement_type == "foldback")
      ) +
      ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
        alpha = svalpha,
        size = svwidth,
        col = "grey30",
        data = bezdf %>% dplyr::filter(rearrangement_type != "foldback")
      )
  }
  
  if (!is.null(genes)) {
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    gene_idx <- get_gene_idx(genes, chr = chrfilt, binsize = binsize)
    npoints <- dim(pl$CNbins)[1]
    gCN <- gCN +
      ggplot2::geom_vline(data = gene_idx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3) +
      ggrepel::geom_text_repel(data = gene_idx, ggplot2::aes(x = idx - npoints * adj, y = maxCN, label = ensembl_gene_symbol), col = "black", alpha = 0.75)
  }

  if (!is.null(annotateregions)) {
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    annotateregions <- dplyr::mutate(annotateregions, start = round(start / binsize) * binsize + 1)
    datidx <- dplyr::inner_join(annotateregions, pl$bins %>% dplyr::select(chr, start, idx)) %>% dplyr::distinct(.)
    gCN <- gCN +
      ggplot2::geom_vline(data = datidx, ggplot2::aes(xintercept = idx), lty = annotateregions_linetype, size = 0.3, alpha = 0.5)
  }

  return(gCN)
}

#' Plot a single cell allele specific copy number profile
#'
#' @param cn Single cell allele specific copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy` or a hscn object.
#' @param cellid Which cell to plot, if no cell is specific will plot the first cell in the dataframe
#' @param chrfilt Vector of chromosomes to plot, if NULL (default) will plot all chromosomes
#' @param pointsize The point size in the plot
#' @param alphaval Alpha value of points
#' @param maxCN The maximum on the y axis, if any points are above this value they will be winsorized rather than removed
#' @param cellidx idx of cell to plot if cellid = NULL
#' @param statecol The colour mapping, default is to map colours to the `state` column
#' @param returnlist Return a list rather than the ggplot object
#' @param raster use ggrastr or not, default = FALSE
#' @param y_axis_trans What transformation to use on the y-axis, default is identity, the other option is "squashy" which uses a tanh transformation
#' @param xaxis_order Default is "genome_position"
#' @param legend.position Where to place the legend, default is "bottom"
#' @param annotateregions Dataframe with chr start and end positions to annotate, will draw a dashed vertical line at this position
#' @param annotateregions_linetype Linetype for region annotation, default = 2 (dashed)
#' @param SV Default is NULL. If a dataframe with structural variant position is passed it will add a track on the top showin rearrangement links
#' @param svalpha the alpha scaling of the SV lines, default = 0.5
#' @param svwidth width of rearrangement connections
#' @param genes vector of genes to annotate, will add a dashed vertical line and label
#' @param homolog Rather than plot the BAF and CN seperately this will plot the 2 homologs on the same track
#' @param plotdata Binary value whether to plot raw data or inferred states in homolog plot
#' @param offest to use when plotting inferred states in homolog plot
#' @param chrstart Start of region (in Mb) when plotting a single chromosome
#' @param chrend End of region (in Mb) when plotting a single chromosome
#' @param shape shape for plotting, default = 16
#' @param positionticks set to TRUE to use position ticks rather than chromosome ticks
#' @param BAFcol state to use to colour BAF track, default = `state_phase`
#' @param my_title string to use for title, if NULL cell_id is shown
#'
#'
#' @return ggplot2 plot
#'
#' @examples
#' 
#' \dontrun{
#' data("haplotypes")
#' data("CNbins")
#' haplotypes <- format_haplotypes_dlp(haplotypes, CNbins)
#' hscn <- callHaplotypeSpecificCN(CNbins, haplotypes, likelihood = "binomial")
#' plotCNprofileBAF(hscn, genes = "MYC")
#' plotCNprofileBAF(hscn, homolog = TRUE, chrfilt = c("1", "8"))
#' }
#'
#' @export
plotCNprofileBAF <- function(cn,
                             cellid = NULL,
                             chrfilt = NULL,
                             pointsize = 1,
                             alphaval = 0.6,
                             maxCN = 10,
                             cellidx = 1,
                             BAFcol = "state_phase",
                             statecol = "state",
                             returnlist = FALSE,
                             raster = FALSE,
                             y_axis_trans = "identity",
                             xaxis_order = "genome_position",
                             legend.position = "bottom",
                             genes = NULL,
                             annotateregions = NULL,
                             annotateregions_linetype = 2,
                             homolog = FALSE,
                             SV = NULL,
                             adj = 0.03,
                             svalpha = 0.5,
                             svwidth = 1.0,
                             plotdata = TRUE,
                             offset = NULL,
                             my_title = NULL,
                             tickwidth = 50,                          
                             chrstart = NULL,
                             chrend = NULL,
                             shape = 16,
                             positionticks = FALSE,
                             ...) {
  if (homolog == TRUE) {
    ghomolog <- plotCNprofileBAFhomolog(cn,
      cellid = cellid,
      chrfilt = chrfilt,
      pointsize = pointsize,
      alphaval = alphaval,
      maxCN = maxCN,
      raster = raster,
      cellidx = cellidx,
      returnlist = returnlist,
      y_axis_trans = y_axis_trans,
      xaxis_order = xaxis_order,
      legend.position = legend.position,
      genes = genes,
      annotateregions = annotateregions,
      annotateregions_linetype = annotateregions_linetype,
      SV = SV,
      adj = adj,
      svalpha = svalpha,
      svwidth = svwidth,
      plotdata = plotdata,
      offset = offset,
      tickwidth = tickwidth,
      chrstart = chrstart,
      chrend = chrend,
      shape = shape,
      positionticks = positionticks,
      ...
    )
    return(ghomolog)
  }

  if (!xaxis_order %in% c("bin", "genome_position")) {
    stop("xaxis_order must be either 'bin' or 'genome_position'")
  }

  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }

  if (y_axis_trans == "squashy") {
    ybreaks <- c(0, 2, 5, 10, maxCN)
  } else {
    ybreaks <- seq(0, maxCN, 2)
  }
  
  if (is.null(my_title)){
    my_title <- cellid
  }
  
  if (length(chrfilt) == 1){
    xlab <- paste0('Chr. ', chrfilt, " (Mb)")
  } else{
    xlab <- 'Chromosome'
  }

  if (is.null(cellid)) {
    cellid <- unique(CNbins$cell_id)[min(cellidx, length(unique(CNbins$cell_id)))]
  }

  if (!"BAF" %in% names(CNbins)) {
    stop("No BAF column in dataframe, first calculate the BAF per bin using combineBAFCN and then callAlleleSpecificCN")
  }

  if (BAFcol == "state_min") {
    BAFcolpal <- scCNminorallele_cols()
  }

  if (BAFcol == "state_phase") {
    BAFcolpal <- scCNphase_cols()
  }

  if (BAFcol == "state_AS") {
    BAFcolpal <- scCNAS_cols()
  }

  statecolpal <- scCNstate_cols()

  message(paste0("Making CN profile and BAF plot for cell - ", cellid))

  if (!is.null(chrfilt)) {
    message(paste0("Filtering for chromosomes: ", paste0(chrfilt, collapse = ",")))
    CNbins <- dplyr::filter(CNbins, chr %in% chrfilt)
  }

  pl <- CNbins %>%
    dplyr::filter(cell_id == cellid) %>%
    plottinglist(., xaxis_order = xaxis_order, maxCN = maxCN, positionticks = positionticks,
                 tickwidth = tickwidth, chrstart = chrstart, chrend = chrend)

  if (raster == TRUE) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package \"ggrastr\" needed for this function to work. Please install it.",
        call. = FALSE
      )
    }
    gBAF <- pl$CNbins %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = BAF)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = BAFcol), size = pointsize, alpha = alphaval, shape = shape) +
      ggplot2::scale_color_manual(
        name = "CN",
        breaks = names(BAFcolpal),
        labels = names(BAFcolpal),
        values = BAFcolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.0)) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab("BAF") +
      ggplot2::ggtitle(my_title) +
      cowplot::theme_cowplot(...) +
      ggplot2::geom_hline(yintercept = 0.5, lty = 2, alpha = 0.5) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha = 1, size = 3, shape = 15)))

    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggrastr::geom_point_rast(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = shape) +
      ggplot2::scale_color_manual(
        name = "Allele Specific CN",
        breaks = names(statecolpal),
        labels = names(statecolpal),
        values = statecolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot(...) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = 6, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3, shape = 15)
      )) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
  } else {
    gBAF <- pl$CNbins %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = BAF)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggplot2::geom_point(ggplot2::aes_string(col = BAFcol), size = pointsize, alpha = alphaval, shape = shape) +
      ggplot2::scale_color_manual(
        name = "CN",
        breaks = names(BAFcolpal),
        labels = names(BAFcolpal),
        values = BAFcolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), limits = c(0, 1.0)) +
      ggplot2::xlab("Chromosome") +
      ggplot2::ylab("BAF") +
      ggplot2::ggtitle(my_title) +
      cowplot::theme_cowplot(...) +
      ggplot2::geom_hline(yintercept = 0.5, lty = 2, alpha = 0.5) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      ) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position) +
      ggplot2::guides(colour = ggplot2::guide_legend(ncol = 5, override.aes = list(alpha = 1, size = 3, shape = 15)))

    gCN <- pl$CNbins %>%
      dplyr::mutate(state = ifelse(state >= 11, "11+", paste0(state))) %>%
      dplyr::mutate(state = factor(paste0(state), levels = c(paste0(seq(0, 10, 1)), "11+"))) %>%
      dplyr::mutate(state_min = paste0(state_min)) %>%
      ggplot2::ggplot(ggplot2::aes(x = idx, y = copy)) +
      ggplot2::geom_vline(xintercept = pl$chrbreaks, col = "grey90", alpha = 0.75) +
      ggplot2::geom_point(ggplot2::aes_string(col = statecol), size = pointsize, alpha = alphaval, shape = shape) +
      ggplot2::scale_color_manual(
        name = "Allele Specific CN",
        breaks = names(statecolpal),
        labels = names(statecolpal),
        values = statecolpal,
        drop = FALSE
      ) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "none"
      ) +
      ggplot2::scale_x_continuous(breaks = pl$chrticks, labels = pl$chrlabels, expand = c(0, 0), limits = c(pl$minidx, pl$maxidx), guide = ggplot2::guide_axis(check.overlap = TRUE)) +
      ggplot2::scale_y_continuous(breaks = ybreaks, limits = c(0, maxCN), trans = y_axis_trans) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab("Copy Number") +
      cowplot::theme_cowplot(...) +
      ggplot2::guides(colour = ggplot2::guide_legend(
        ncol = 6, byrow = TRUE,
        override.aes = list(alpha = 1, size = 3, shape = 15)
      )) +
      ggplot2::theme(legend.title = ggplot2::element_blank(), legend.position = legend.position)
  }

  if (!is.null(SV)) {
    svpl <- plottinglistSV(SV, chrfilt = chrfilt)
    bezdf <- get_bezier_df(svpl, pl, maxCN)
    bezdf <- bezdf %>%
      dplyr::filter((position_1 != position_2) | rearrangement_type == "foldback")
    gCN <- gCN +
      ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
        alpha = 0.8,
        size = svwidth,
        col = as.vector(SV_colors["Foldback"]),
        data = bezdf %>% dplyr::filter(rearrangement_type == "foldback")
      ) +
      ggforce::geom_bezier(ggplot2::aes(x = idx, y = copy, group = id),
        alpha = svalpha,
        col = "grey30",
        size = svwidth,
        data = bezdf %>% dplyr::filter(rearrangement_type != "foldback")
      )
  }

  if (!is.null(genes)) {
    binsize <- pl$CNbins$end[1] - pl$CNbins$start[1] + 1
    gene_idx <- get_gene_idx(genes, chr = chrfilt, binsize = binsize)
    npoints <- dim(pl$CNbins)[1]
    gBAF <- gBAF +
      ggplot2::geom_vline(data = gene_idx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3) +
      ggrepel::geom_text_repel(data = gene_idx, ggplot2::aes(x = idx - npoints * adj, y = 1.0, label = ensembl_gene_symbol), col = "black", alpha = 0.75)
    gCN <- gCN +
      ggplot2::geom_vline(data = gene_idx, ggplot2::aes(xintercept = idx), lty = 2, size = 0.3)
  }

  if (!is.null(annotateregions)) {
    datidx <- dplyr::inner_join(annotateregions, pl$bins %>% dplyr::select(chr, start, idx)) %>% dplyr::distinct(.)
    gBAF <- gBAF +
      ggplot2::geom_vline(data = datidx, ggplot2::aes(xintercept = idx), lty = annotateregions_linetype, size = 0.3, alpha = 0.5)
    gCN <- gCN +
      ggplot2::geom_vline(data = datidx, ggplot2::aes(xintercept = idx), lty = annotateregions_linetype, size = 0.3, alpha = 0.5)
  }

  g <- cowplot::plot_grid(gBAF, gCN, align = "v", ncol = 1, rel_heights = c(1, 1.2))

  if (returnlist == TRUE) {
    p <- list(CN = gCN, BAF = gBAF, both = g, plist = pl)
  } else {
    p <- g
  }

  return(p)
}


#' @export
plotCNBAF <- function(cn, nfilt = 10^5, plottitle = "5Mb", pointsize = 0.1, shape = 16, ...) {
  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }
  maj <- seq(1, 11, 1)
  min <- seq(0, 11, 1)
  ASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS = paste0(state - min, "|", min)) %>%
    dplyr::mutate(A = state - min, B = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(A >= 0, B >= 0)

  CNbins <- CNbins %>%
    dplyr::filter(state < 10, copy < 10)

  if (nfilt < length(CNbins$chr)) {
    CNbins <- CNbins %>%
      dplyr::sample_n(nfilt)
  }

  plottitle <- paste0((CNbins$end[1] - CNbins$start[1] + 1) / 1e6, "Mb bins")

  g <- CNbins %>%
    ggplot2::ggplot(ggplot2::aes(y = BAF, x = copy, col = paste0("CN", state))) +
    ggplot2::geom_point(size = pointsize, alpha = 0.2, shape = shape) +
    ggplot2::xlab("Corrected read counts") +
    cowplot::theme_cowplot(...) +
    ggplot2::ggtitle(plottitle) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 5))) +
    ggplot2::geom_text(data = ASstates, ggplot2::aes(x = state, y = cBAF, label = state_AS), col = "black") +
    ggplot2::scale_color_manual(
      name = "Copy number \n state",
      breaks = paste0("CN", seq(0, max(CNbins$state, na.rm = TRUE), 1)),
      labels = seq(0, max(CNbins$state, na.rm = TRUE), 1),
      values = scCN_cols(paste0("CN", seq(0, max(CNbins$state, na.rm = TRUE), 1)))
    ) +
    ggplot2::scale_x_continuous(breaks = seq(0, 10, 1), labels = seq(0, 10, 1), limits = c(0, 10))
  return(g)
}

#' Plot BAF distributions per allele specific state
#'
#' @param cn Either a hscn object or a single cell allele specific copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`, `BAF`
#' @param minpts minimum number of points to include default = 250
#' @param minfrac states that are present below this fraction will be removed, default = 0.01
#' @param maxstate States with total copy number > maxstate will be removed, default = 10
#' @param dens_adjust density adjustment factor in the violin plots
#' @param mincounts filter out bins < mincounts from plotting, default = 6
#'
#' @examples
#' \dontrun{
#' data("haplotypes")
#' data("CNbins")
#' haplotypes <- format_haplotypes_dlp(haplotypes, CNbins)
#' hscn <- callHaplotypeSpecificCN(CNbins, haplotypes, likelihood = "binomial")
#' plotBAFperstate(hscn, genes = "MYC")
#' }
#'
#' @export
plotBAFperstate <- function(cn, minpts = 250, minfrac = 0.01, maxstate = 10, dens_adjust = 2.0, mincounts = 6) {
  if (is.hscn(cn) | is.ascn(cn)) {
    alleleCN <- cn$data
  } else {
    alleleCN <- cn
  }

  maj <- seq(0, max(alleleCN$state), 1)
  min <- seq(0, max(alleleCN$state), 1)

  allASstates <- expand.grid(state = maj, min = min) %>%
    dplyr::mutate(cBAF = min / state) %>%
    dplyr::mutate(state_AS_phased = paste0(state - min, "|", min)) %>%
    dplyr::mutate(A = state - min, B = min) %>%
    dplyr::select(-min) %>%
    dplyr::filter(A >= 0, B >= 0) %>%
    dplyr::filter(state <= maxstate)
  allASstates$cBAF[is.nan(allASstates$cBAF)] <- 0.0

  forplot <- alleleCN %>%
    dplyr::filter(totalcounts > mincounts) %>% 
    dplyr::group_by(state_AS_phased) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(f = n / dplyr::n()) %>%
    dplyr::filter(state <= maxstate, n > minpts, f > minfrac)
  allASstates <- allASstates %>%
    dplyr::filter(state_AS_phased %in% unique(forplot$state_AS_phased))

  text_fraction <- dplyr::distinct(forplot, state_AS_phased, f) %>%
    dplyr::mutate(pct = paste0(100 * round(f, 2), "%"), y = 0.95)
  
  maxCN <- min(11, max(alleleCN$state, na.rm = TRUE))
  minCN <- min(alleleCN$state, na.rm = TRUE)

  g <- forplot %>%
    dplyr::mutate(cncol = paste0("CN", state)) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = forcats::fct_reorder(state_AS_phased, state),
      y = BAF
    )) +
    ggplot2::geom_violin(scale = "width", col = "white", ggplot2::aes(fill = cncol), adjust = dens_adjust) +
    ggplot2::geom_boxplot(width = 0.1, outlier.shape = NA, col = "white", ggplot2::aes(fill = cncol)) +
    ggplot2::scale_fill_manual(
      name = "Copy number \n state",
      breaks = paste0("CN", seq(minCN, maxCN, 1)),
      labels = seq(minCN, maxCN, 1),
      values = scCN_cols(paste0("CN", seq(minCN, maxCN, 1)))
    ) +
    cowplot::theme_cowplot() +
    ggplot2::geom_crossbar(
      data = allASstates, ggplot2::aes(y = cBAF, ymin = cBAF, ymax = cBAF),
      alpha = 0.2, size = 0.2
    ) +
    ggplot2::geom_text(data = text_fraction, ggplot2::aes(x = state_AS_phased, y = y, label = pct)) +
    ggplot2::xlab("") +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(g)
}


plot_density_histogram <- function(dat, mystate, rho, nbins = 30, frac = "NA", loherror = 0.02) {
  dat <- dat %>%
    dplyr::filter(state_AS_phased == mystate, totalcounts > 9)

  expBAF <- strsplit(mystate, "\\|")[[1]]
  expBAF <- as.numeric(expBAF[2]) / (as.numeric(expBAF[1]) + as.numeric(expBAF[2]))
  if (expBAF > 0.5) {
    expBAF <- expBAF - loherror
  } else if (expBAF < 0.5) {
    expBAF <- expBAF + loherror
  } else {
    expBAF <- expBAF
  }

  BAF_bb <- VGAM::rbetabinom(length(dat$totalcounts),
    size = dat$totalcounts,
    prob = expBAF,
    rho = rho
  ) / dat$totalcounts
  dffit_bb <- data.frame(
    BAF = BAF_bb,
    type = paste("Beta Binomial (rho = ", round(rho, 4), ")", sep = "")
  )

  BAF_b <- VGAM::rbetabinom(length(dat$totalcounts),
    size = dat$totalcounts,
    prob = expBAF,
    rho = 0.0
  ) / dat$totalcounts
  dffit_b <- data.frame(
    BAF = BAF_b,
    type = paste("Binomial", sep = "")
  )

  g <- ggplot2::ggplot(dat, ggplot2::aes(x = BAF)) +
    ggplot2::geom_histogram(ggplot2::aes(y = (..density..)), bins = nbins, fill = "azure4", alpha = 0.8) +
    ggplot2::geom_line(data = dffit_bb, stat = "density", ggplot2::aes(BAF, col = type), size = 0.5, adjust = 5) +
    ggplot2::geom_line(data = dffit_b, stat = "density", ggplot2::aes(BAF, col = type), size = 0.5, adjust = 5, linetype = 2) +
    ggplot2::scale_colour_manual(values = c(ggplot2::alpha("deepskyblue4", 0.6), ggplot2::alpha("firebrick4", 0.6))) +
    ggplot2::theme_bw() +
    ggplot2::xlab("BAF") +
    ggplot2::ylab("Density") +
    ggplot2::xlim(c(0.0, 1.0)) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    cowplot::theme_cowplot() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::ggtitle(paste0(mystate, " (", round(frac, 3) * 100, "%)"))

  return(g)
}

#' @export
plotBBfit <- function(hscn, nbins = 30, minfrac = 0.01) {
  x <- table(hscn$data$state_AS_phased)
  x <- x / sum(x)
  x <- x[which(x > minfrac)]

  mydat <- dplyr::filter(hscn$data, A != 0, B != 0, totalcounts > 9)

  x <- x[names(x) %in% unique(mydat$state_AS_phased)]

  gplots <- list()
  j <- 1
  for (i in names(x)) {
    gplots[[j]] <- plot_density_histogram(hscn$data,
      mystate = i,
      rho = hscn$likelihood$rho,
      nbins = nbins,
      frac = x[[i]],
      loherror = hscn$loherror
    )
    j <- j + 1
  }

  legend <- cowplot::get_legend(
    # create some space to the left of the legend
    gplots[[j - 1]] + ggplot2::theme(legend.box.margin = ggplot2::margin(0, 0, 0, 12)) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) + ggplot2::theme(legend.position = "right")
  )

  gplots[[j]] <- legend

  cowplot::plot_grid(plotlist = gplots, ncol = 3)
}

#' @export
plot_variance_state <- function(hscn, by_allele_specific_state = FALSE) {
  if (is.hscn(hscn) | is.ascn(hscn)) {
    dat <- hscn$data
  } else {
    dat <- hscn
  }

  if (by_allele_specific_state == TRUE) {
    plot_var <- dat %>%
      dplyr::group_by(state, state_AS_phased) %>%
      dplyr::filter(B > 0 & A > 0) %>%
      dplyr::summarize(mBAF = median(BAF), varBAF = var(BAF)) %>%
      ggplot2::ggplot(ggplot2::aes(x = state, y = varBAF, col = factor(paste0("CN", state), levels = paste0("CN", seq(0, 11, 1))))) +
      ggplot2::geom_text(ggplot2::aes(label = state_AS_phased)) +
      cowplot::theme_cowplot() +
      ggplot2::xlab("State") +
      ggplot2::ylab("BAF variance") +
      ggplot2::scale_color_manual(
        name = "Copy number \n state",
        breaks = paste0("CN", seq(0, 11, 1)),
        labels = paste0("CN", seq(0, 11, 1)),
        values = scCN_cols(paste0("CN", seq(0, 11, 1))),
        drop = FALSE
      ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(
        labels = seq(1, max(dat$state), 1),
        breaks = seq(1, max(dat$state), 1)
      )
  } else {
    plot_var <- dat %>%
      dplyr::filter(B > 0 & A > 0) %>%
      dplyr::group_by(state, state_AS_phased) %>%
      dplyr::summarize(mBAF = median(BAF), varBAF = var(BAF)) %>%
      dplyr::group_by(state) %>%
      dplyr::summarise(mBAF = mean(mBAF), ymin = min(varBAF), ymax = max(varBAF), varBAF = mean(varBAF)) %>%
      ggplot2::ggplot(ggplot2::aes(
        x = state, y = varBAF,
        ymin = ymin, ymax = ymax,
        col = factor(paste0("CN", state), levels = paste0("CN", seq(0, 11, 1)))
      )) +
      ggplot2::geom_pointrange() +
      cowplot::theme_cowplot() +
      ggplot2::xlab("State") +
      ggplot2::ylab("BAF variance") +
      ggplot2::scale_color_manual(
        name = "Copy number \n state",
        breaks = paste0("CN", seq(0, 11, 1)),
        labels = paste0("CN", seq(0, 11, 1)),
        values = scCN_cols(paste0("CN", seq(0, 11, 1))),
        drop = FALSE
      ) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_x_continuous(
        labels = seq(1, max(dat$state), 1),
        breaks = seq(1, max(dat$state), 1)
      )
  }

  return(plot_var)
}

#' Plot clusters used for phasing
#'
#' Plots the consensus copy number and BAF profiles of clusters used for phasing
#'
#' @param hscn An object of class \code{hscn} generated by \code{\link{hscn}}
#'
#' @return A \code{ggplot2} object
#'
#'
#' @export
plot_clusters_used_for_phasing <- function(hscn){
  myplots <- list()
  for (mychr in gtools::mixedsort(names(hscn$phasing))){
    cells <- hscn$phasing[[mychr]]
    myplots[[mychr]] <- hscn$data %>% 
      dplyr::filter(cell_id %in% cells) %>% 
      consensuscopynumber(.) %>% 
      dplyr::mutate(cell_id = paste0("chr ", mychr)) %>% 
      plotCNprofileBAF(., chrfilt = mychr, legend.position = "none")
  }
  
  g <- cowplot::plot_grid(plotlist = myplots, ncol = 6)
  return(g)
}