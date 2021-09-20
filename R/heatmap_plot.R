
cn_colours <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
)

cn_colours_loh <- scCNAS_colors
cn_colours_minorallele <- scCNminorallele_colors
cn_colours_phase <- scCNphase_colors
cn_colours_bafstate <- scBAFstate_colors

snv_colours <- structure(
  names = c(0, 1),
  c("#7EA5EA", "#9A2E1C")
)

clone_palette_20 <- c(
  "#be5f72", "#d74058", "#dc4229", "#a6552c", "#df956f", "#e47a33",
  "#d49f34", "#836e2c", "#b2ad5a", "#92b539", "#4c7d38", "#4dc041",
  "#5dba7f", "#47b8c3", "#6280ca", "#7b57db", "#ce8bd1", "#934f94",
  "#cb48cb", "#d74391"
)
clone_none_black <- "#1B1B1B"

calc_state_mode <- function(states) {
  state_levels <- unique(states)
  state_mode <- state_levels[
    which.max(tabulate(match(states, state_levels)))
  ]
  if (!is.finite(state_mode)) {
    state_mode <- 2
  }
  return(state_mode)
}

normalize_cell_ploidy <- function(copynumber) {
  cell_ids <- colnames(copynumber)
  cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

  for (cell_id in cell_ids) {
    state_mode <- calc_state_mode(copynumber[[cell_id]])
    copynumber[[cell_id]] <- as.integer(ceiling(
      copynumber[[cell_id]] / (state_mode / 2)
    ))
    copynumber[[cell_id]][copynumber[[cell_id]] > 11] <- 11
  }
  return(copynumber)
}

createSVmatforhmap <- function(x, cnmat) {
  options(scipen = 999)

  breakends <- dplyr::bind_rows(
    dplyr::select(x, chromosome_1, position_1, rearrangement_type, read_count) %>% dplyr::rename(chromosome = chromosome_1, position = position_1),
    dplyr::select(x, chromosome_2, position_2, rearrangement_type, read_count) %>% dplyr::rename(chromosome = chromosome_2, position = position_2)
  )

  breakends <- breakends %>%
    dplyr::mutate(position = 0.5e6 * floor(position / 0.5e6) + 1) %>%
    dplyr::mutate(loci = paste0(chromosome, ":", position, ":", position + 0.5e6 - 1)) %>%
    dplyr::arrange(desc(read_count)) %>%
    dplyr::group_by(loci) %>%
    dplyr::filter(row_number() == 1) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(loci, rearrangement_type)

  breakends_ <- data.frame(loci = names(cnmat)) %>%
    dplyr::left_join(breakends) %>%
    dplyr::mutate(y = ifelse(is.na(rearrangement_type), 0, 1)) %>%
    dplyr::mutate(col = case_when(
      is.na(rearrangement_type) ~ NA_character_,
      rearrangement_type == "inversion" ~ SV_colors[["Inversion"]],
      rearrangement_type == "foldback" ~ SV_colors[["Foldback"]],
      rearrangement_type == "unbalanced" ~ SV_colors[["Unbalanced"]],
      rearrangement_type == "duplication" ~ SV_colors[["Duplication"]],
      rearrangement_type == "balanced" ~ SV_colors[["Balanced"]],
      rearrangement_type == "deletion" ~ SV_colors[["Deletion"]]
    ))

  return(breakends_)
}

format_tree <- function(tree, brlen) {
  locus_tips <- grep("locus", tree$tip.label, value = TRUE)
  tree <- ape::drop.tip(tree, locus_tips)

  if (!is.null(brlen)) {
    tree <- ape::compute.brlen(tree, brlen)
  }

  tree$tip.label <- gsub("cell_", "", tree$tip.label)

  return(tree)
}

cell_order_from_tree <- function(tree, clones, cells = NULL) {
  # tree <- ape::ladderize(tree, right = TRUE)
  cellorder <- tree$tip.label[tree$edge[, 2]]
  clones <- as.data.frame(clones)
  row.names(clones) <- clones$cell_id
  clones <- clones[cellorder, ]

  if (!is.null(cells)) {
    clones <- clones[clones$cell_id %in% cells, ]
  }
  clones <- clones[!is.na(clones$cell_id), ]
  return(clones)
}

get_clone_members <- function(clones) {
  clone_members <- list()
  for (c in unique(clones$clone_id)) {
    if (c != "None") {
      clone_members[[c]] <- clones[clones$clone_id == c, "cell_id"]
    }
  }
  return(clone_members)
}

make_clone_palette <- function(levels) {
  if (length(levels) <= 8) {
    pal <- RColorBrewer::brewer.pal(max(length(levels), 3), "Paired")
  } else {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), "Paired"))(length(levels))
  }
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

make_tree_ggplot <- function(tree, clones, clone_pal = NULL) {
  if (!is.null(clones)) {
    clone_members <- get_clone_members(clones)
    tree <- ggtree::groupOTU(tree, clone_members)
    clone_levels <- gtools::mixedsort(unique(clones$clone_id))
    if (is.null(clone_pal)) {
      clone_pal <- make_clone_palette(clone_levels)
    }
    clone_pal[["0"]] <- clone_none_black
    tree_aes <- ggplot2::aes(x, y, colour = group)
  } else {
    tree_aes <- ggplot2::aes(x, y)
  }

  p <- ggplot2::ggplot(tree, tree_aes) +
    ggtree::geom_tree(size = 0.25) +
    ggplot2::coord_cartesian(expand = FALSE) +
    ggplot2::ylim(0.5, length(tree$tip.label) + 0.5) +
    ggplot2::theme_void()

  if (!is.null(clones)) {
    p <- p + ggplot2::scale_colour_manual(values = clone_pal)
  }

  return(p)
}

make_discrete_palette <- function(pal_name, levels) {
  if (length(levels) > 8) {
    pal_name <- "Set3"
  }
  pal <- colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), pal_name))(length(levels))
  # pal <- RColorBrewer::brewer.pal(max(length(levels), 3), pal_name)
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

format_copynumber_values <- function(copynumber, plotcol = "state") {
  # copynumber[copynumber > 11] <- 11

  if (plotcol %in% c("BAF", "copy")) {
    for (col in colnames(copynumber)) {
      values <- copynumber[, col]
      copynumber[, col] <- values
    }
  } else {
    for (col in colnames(copynumber)) {
      values <- as.character(copynumber[, col])
      values[values == "11"] <- "11+"
      copynumber[, col] <- values
    }
  }
  return(copynumber)
}

space_copynumber_columns <- function(copynumber, spacer_cols) {
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  spacer <- as.data.frame(matrix(
    data = NA, nrow = nrow(copynumber), ncol = spacer_cols
  ))
  chrom_copynumber_dfs <- list()
  for (chrom in gtools::mixedsort(unique(chroms))) {
    chrom_copynumber <- copynumber[, chroms == chrom, drop = FALSE]
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(chrom_copynumber))
    chrom_copynumber_dfs <- c(chrom_copynumber_dfs, list(spacer))
  }
  chrom_copynumber_dfs[length(chrom_copynumber_dfs)] <- NULL
  copynumber <- do.call(cbind, chrom_copynumber_dfs)

  return(copynumber)
}

get_ordered_cell_ids <- function(tree_plot_dat) {
  return(rev(dplyr::arrange(tree_plot_dat[tree_plot_dat$isTip, ], y)$label))
}

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE) {
  do.call(order, c(
    lapply(list(...), function(l) {
      if (is.character(l)) {
        factor(l, levels = gtools::mixedsort(unique(l)))
      } else {
        l
      }
    }),
    list(na.last = na.last, decreasing = decreasing)
  ))
}

format_copynumber <- function(copynumber,
                              ordered_cell_ids,
                              plotcol = "state",
                              spacer_cols = 20) {
  if (!("chr" %in% colnames(copynumber))) {
    message("No chr column")
    loci <- sapply(rownames(copynumber), strsplit, "_")
    copynumber$chr <- unname(sapply(loci, "[[", 1))
    copynumber$start <- as.numeric(unname(sapply(loci, "[[", 2)))
    copynumber$end <- as.numeric(unname(sapply(loci, "[[", 3)))
    copynumber$width <- (copynumber$end - copynumber$start + 1)
  }
  copynumber$chr <- gsub("chr", "", copynumber$chr)
  copynumber <- copynumber[multi.mixedorder(copynumber$chr, copynumber$start), ]

  rownames(copynumber) <- paste0(
    copynumber$chr, ":", copynumber$start, ":", copynumber$end
  )
  copynumber <- subset(copynumber, select = -c(chr, start, end, width))
  copynumber <- as.data.frame(t(copynumber))

  copynumber <- copynumber[ordered_cell_ids, ]

  copynumber <- format_copynumber_values(copynumber, plotcol = plotcol)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)

  return(copynumber)
}

format_clones <- function(clones, ordered_cell_ids) {
  clonesdf <- dplyr::full_join(clones, data.frame(cell_id = ordered_cell_ids), by = "cell_id")
  clonesdf[is.na(clonesdf$clone_id), "clone_id"] <- "None"
  clone_counts <- clones %>%
    dplyr::group_by(clone_id) %>%
    dplyr::summarise(count = dplyr::n())

  clonesdf <- dplyr::left_join(clonesdf, clone_counts, by = "clone_id") %>%
    dplyr::mutate(clone_label = paste0(clone_id, " (", count, ")")) %>%
    dplyr::select(-count) %>%
    as.data.frame()

  rownames(clonesdf) <- clonesdf$cell_id
  clonesdf <- clonesdf[ordered_cell_ids, ]
  return(clonesdf)
}

make_corrupt_tree_heatmap <- function(tree_ggplot, ...) {
  tree_annot_func <- ComplexHeatmap::AnnotationFunction(
    fun = function(index) {
      pushViewport(viewport(height = 1))
      grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[5]])
      popViewport()
    },
    var_import = list(tree_ggplot = tree_ggplot),
    width = grid::unit(4, "cm"),
    which = "row"
  )
  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree = tree_annot_func, which = "row", show_annotation_name = FALSE
  )

  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- ComplexHeatmap::Heatmap(matrix(ncol = 0, nrow = n_cells), left_annotation = tree_annot, ...)

  return(tree_hm)
}

find_largest_contiguous_group <- function(x) {
  starts <- c(1, which(diff(x) != 1 & diff(x) != 0) + 1)
  ends <- c(starts[-1] - 1, length(x))
  largest <- which.max(ends - starts + 1)
  return(x[starts[largest]:ends[largest]])
}

get_clone_label_pos <- function(clones) {
  clone_label_pos <- list()
  for (clone in unique(clones$clone_id)) {
    if (!grepl("None", clone)) {
      clone_idx <- which(clones$clone_id == clone)
      clone_idx <- find_largest_contiguous_group(clone_idx)
      clone_label_pos[[as.character(clone)]] <-
        as.integer(round(mean(clone_idx)))
    }
  }
  return(clone_label_pos)
}

get_label <- function(cell_id, idx, str_to_remove){
  totndash <- stringr::str_count(cell_id, "-")
  ndash <- totndash - 2
  #remove last 2 dashes (-R*-C*)
  lab <- sub(paste0("^(([^-]*-){", ndash, "}[^-]*).*"), "\\1", cell_id)
  if (idx == 2){
    #get the library id
    lab <- strsplit(lab, "-")[[1]][ndash + 1]
  } else if (idx == 1){
    #get the sample id (complicated because some sample IDs have -'s)
    lab <- sub(paste0("^(([^-]*-){", ndash - 1, "}[^-]*).*"), "\\1", lab)
  }
  
  if (!is.null(str_to_remove)){
    lab <- stringr::str_remove(lab, str_to_remove)
  }
  
  return(lab)
}

get_library_labels <- function(cell_ids, idx = 1, str_to_remove = NULL) {
  labels <- sapply(cell_ids, function(x) {
    return(get_label(x, idx, str_to_remove))
  })
  return(labels)
}

make_left_annot <- function(copynumber,
                            clones,
                            library_mapping = NULL,
                            show_library_label = TRUE,
                            show_clone_label = TRUE,
                            clone_pal = NULL,
                            idx = 1,
                            show_legend = TRUE,
                            str_to_remove = NULL) {
  annot_colours <- list()

  if (show_clone_label == FALSE & show_library_label == FALSE) {
    return(NULL)
  }

  library_labels <- get_library_labels(rownames(copynumber), idx = idx, str_to_remove = str_to_remove)
  if (!is.null(library_mapping)) {
    library_labels <- unlist(library_mapping[library_labels])
    if (!all(library_labels %in% names(library_mapping)) == FALSE) {
      warning("Not all library ids present in library to name mapping, using library IDs as annotations...")
      library_labels <- get_library_labels(rownames(copynumber))
    }
  }
  library_levels <- gtools::mixedsort(unique(library_labels))
  annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  annot_colours$Sample <- annot_colours$Sample[!is.na(annot_colours$Sample)]

  library_legend_rows <- 3

  if (!is.null(clones)) {
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])
    if (is.null(clone_pal)) {
      clone_pal <- make_clone_palette(clone_levels)
    }

    if (length(clone_level_none > 0)) {
      clone_pal[[clone_level_none]] <- clone_none_black
    }
    annot_colours$Cluster <- clone_pal

    clone_label_generator <- function(index) {
      clone_label_pos <- get_clone_label_pos(clones)
      y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
      grid::grid.text(
        names(clone_label_pos), 0.5, y_pos,
        just = c("centre", "centre")
      )
    }

    clone_legend_rows <- 3
    if (length(clone_levels) > 3) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 2))
    }

    if (show_library_label == TRUE & show_clone_label == TRUE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Cluster = clones$clone_label, clone_label = clone_label_generator,
        Sample = library_labels,
        col = annot_colours, show_annotation_name = c(TRUE, FALSE, TRUE),
        which = "row", annotation_width = grid::unit(rep(0.4, 3), "cm"),
        annotation_legend_param = list(
          Cluster = list(nrow = clone_legend_rows, direction = "horizontal"),
          Sample = list(nrow = library_legend_rows, direction = "horizontal")
        ),
        show_legend = show_legend
      )
    } else if (show_library_label == FALSE & show_clone_label == TRUE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Cluster = clones$clone_label, clone_label = clone_label_generator,
        col = annot_colours, show_annotation_name = c(TRUE, FALSE),
        which = "row", annotation_width = grid::unit(rep(0.4, 2), "cm"),
        annotation_legend_param = list(
          Cluster = list(nrow = clone_legend_rows)
        ),
        show_legend = show_legend
      )
    } else if (show_library_label == TRUE & show_clone_label == FALSE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Sample = library_labels, col = annot_colours,
        which = "row", simple_anno_size = grid::unit(0.4, "cm"),
        annotation_legend_param = list(
          Sample = list(nrow = library_legend_rows)
        ),
        show_legend = show_legend
      )
    }
  } else {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Sample = library_labels, col = annot_colours,
      which = "row", simple_anno_size = grid::unit(0.4, "cm"),
      annotation_legend_param = list(
        Sample = list(nrow = library_legend_rows)
      ),
      show_legend = show_legend
    )
  }

  return(left_annot)
}

make_top_annotsnv <- function(mutgroups) {

}

get_genomecoords_label_pos <- function(copynumber, nticks = 3) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  binwidth <- as.numeric(strsplit(colnames(copynumber)[1], ":")[[1]][3]) -
    as.numeric(strsplit(colnames(copynumber)[1], ":")[[1]][2]) + 1
  chromfreq <- table(chroms)
  uniq_chroms <- names(chromfreq)[chromfreq > 1]
  uniq_chroms <- uniq_chroms[stringr::str_detect(uniq_chroms, "V", negate = TRUE)]
  for (chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    max_idx <- round(0.9 * max(chrom_idx))
    nwidth <- 10 * round((max(chrom_idx) / nticks) / 10)
    mypos <- nwidth
    for (ticks in 1:nticks){
      if (mypos > max_idx){
        next
      }
      lab <- paste0(mypos * binwidth / 1e6, "Mb")
      chrom_label_pos[[lab]] <- mypos
      mypos <- mypos + nwidth
    }
  }
  return(chrom_label_pos)
}

get_chrom_label_pos <- function(copynumber, nticks = 3) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  chromfreq <- table(chroms)
  uniq_chroms <- names(chromfreq)[chromfreq > 1]
  uniq_chroms <- uniq_chroms[stringr::str_detect(uniq_chroms, "V", negate = TRUE)]

  if (length(uniq_chroms) == 1){
   chrom_label_pos <- get_genomecoords_label_pos(copynumber, nticks = nticks)
   return(chrom_label_pos)
  }

  for (chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(mean(chrom_idx)))
  }
  return(chrom_label_pos)
}

make_bottom_annot <- function(copynumber,
                              chrlabels = TRUE,
                              filterlabels = NULL,
                              nticks = 3,
                              annotation_height = NULL, 
                              annofontsize = 14,
                              linkheight = 1) {
  if (chrlabels[1] == FALSE) {
    return(NULL)
  } else if (chrlabels[1] == TRUE) {
    chrom_label_pos <- get_chrom_label_pos(copynumber, nticks = nticks)
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels = ComplexHeatmap::anno_mark(
      at = as.vector(unlist(chrom_label_pos)),
      labels = names(chrom_label_pos),
      link_height = grid::unit(linkheight, "mm"),
      labels_gp = grid::gpar(fontsize = annofontsize),
      side = "bottom",
      padding = grid::unit(1.25, "mm"), extend = 0.01, labels_rot = 0
    ), show_annotation_name = FALSE,
    annotation_height = annotation_height)
  } else {
    chrom_label_pos <- get_chrom_label_pos(copynumber)
    chrom_label_pos <- chrom_label_pos[chrlabels]
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels = ComplexHeatmap::anno_mark(
      at = as.vector(unlist(chrom_label_pos)),
      labels = names(chrom_label_pos),
      side = "bottom",
      labels_gp = grid::gpar(fontsize = annofontsize),
      padding = grid::unit(1.25, "mm"), extend = 0.01, labels_rot = 0
    ), show_annotation_name = FALSE,
    annotation_height = annotation_height)
  }
  return(bottom_annot)
}

make_top_annotation_gain <- function(copynumber,
                                     plotcol = "state",
                                     plotfrequency = FALSE,
                                     cutoff = NULL,
                                     maxf = NULL,
                                     SV = NULL) {
  ncells <- nrow(copynumber)

  if ((plotcol == "state" | plotcol == "copy") & plotfrequency == TRUE) {
    f1 <- colSums(copynumber > cutoff, na.rm = TRUE) / ncells
    f2 <- -colSums(copynumber < cutoff, na.rm = TRUE) / ncells
    if (is.null(maxf)) {
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01) {
        maxf <- 0.01
      }
    }
    ha2 <- ComplexHeatmap::columnAnnotation(
      dist2 = ComplexHeatmap::anno_barplot(
        f1,
        bar_width = 1,
        gp = grid::gpar(col = "#E34A33", fill = "#E34A33"),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        f2,
        bar_width = 1,
        gp = grid::gpar(col = "#3182BD", fill = "#3182BD"),
        axis_param = list(
          at = c(0.0, -round(maxf / 2, 2), -maxf),
          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm")
    )
  } else if (plotcol == "state_phase" & plotfrequency == TRUE) {
    f1a <- colSums(apply(copynumber, 2, function(x) grepl("A-Gained", x))) / ncells
    f1b <- colSums(apply(copynumber, 2, function(x) grepl("A-Hom", x))) / ncells
    f2a <- -colSums(apply(copynumber, 2, function(x) grepl("B-Gained", x))) / ncells
    f2b <- -colSums(apply(copynumber, 2, function(x) grepl("B-Hom", x))) / ncells
    f1 <- f1a + f1b
    f2 <- f2a + f2b
    if (is.null(maxf)) {
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01) {
        maxf <- 0.01
      }
    }

    ha2 <- ComplexHeatmap::columnAnnotation(
      dist2 = ComplexHeatmap::anno_barplot(
        matrix(data = c(f1a, f1b), ncol = 2),
        bar_width = 1,
        gp = grid::gpar(
          col = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"]),
          fill = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"])
        ),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        matrix(data = c(f2a, f2b), ncol = 2),
        bar_width = 1,
        gp = grid::gpar(
          col = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"]),
          fill = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"])
        ),
        axis_param = list(
          at = c(0, -round(maxf / 2, 2), -maxf),
          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm")
    )
  }
  else if ((plotcol == "state_BAF" | plotcol == "BAF") & plotfrequency == TRUE) {
    f1 <- colSums(copynumber < 0.5, na.rm = TRUE) / ncells
    f2 <- -colSums(copynumber > 0.5, na.rm = TRUE) / ncells
    if (is.null(maxf)) {
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01) {
        maxf <- 0.01
      }
    }
    ha2 <- ComplexHeatmap::columnAnnotation(
      dist2 = ComplexHeatmap::anno_barplot(
        f1,
        bar_width = 1,
        gp = grid::gpar(col = scCNphase_colors["A-Hom"], fill = scCNphase_colors["A-Hom"]),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        f2,
        bar_width = 1,
        gp = grid::gpar(col = scCNphase_colors["B-Hom"], fill = scCNphase_colors["B-Hom"]),
        axis_param = list(
          at = c(0.0, -round(maxf / 2, 2), -maxf),
          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm")
    )
  }
  else {
    ha2 <- NULL
  }

  if (!is.null(SV)) {
    breakends <- createSVmatforhmap(SV, copynumber)
    annotationbreaks <- sort(unique(breakends$rearrangement_type))
    annotationbreaks <- annotationbreaks[!is.na(annotationbreaks)]
    annotationlabels <- unlist(lapply(annotationbreaks, CapStr))
    ha2 <- ComplexHeatmap::HeatmapAnnotation(
      SV = ComplexHeatmap::anno_barplot(breakends$y,
        gp = grid::gpar(col = breakends$col, fill = breakends$col),
        ylim = c(0, 1),
        axis = FALSE,
        # pch = 25,
        border = FALSE
      ),
      which = "column",
      show_annotation_name = TRUE,
      height = grid::unit(0.7, "cm")
    )
  }

  return(ha2)
}

make_copynumber_heatmap <- function(copynumber,
                                    clones,
                                    colvals = cn_colours,
                                    legendname = "Copy Number",
                                    library_mapping = NULL,
                                    clone_pal = NULL,
                                    sample_label_idx = 1,
                                    cutoff = NULL,
                                    maxf = 1.0,
                                    plotcol = "state",
                                    plotfrequency = FALSE,
                                    show_legend = TRUE,
                                    show_library_label = TRUE,
                                    show_clone_label = TRUE,
                                    chrlabels = TRUE,
                                    SV = NULL,
                                    nticks = 4,
                                    annotation_height = NULL, 
                                    annofontsize = 14,
                                    na_col = "white",
                                    linkheight = 5,
                                    str_to_remove = NULL,
                                    ...) {
  
  if (class(colvals) == "function"){
    leg_params <- list(nrow = 3,
                       direction = "vertical")
  } else {
    leg_params <- list(nrow = 3,
                       direction = "vertical",
                       at = names(colvals))
  }
  
  copynumber_hm <- ComplexHeatmap::Heatmap(
    name = legendname,
    as.matrix(copynumber),
    col = colvals,
    na_col = na_col,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    bottom_annotation = make_bottom_annot(copynumber, chrlabels = chrlabels, nticks = nticks, annotation_height = annotation_height, annofontsize = annofontsize, linkheight = linkheight),
    left_annotation = make_left_annot(copynumber, clones,
      library_mapping = library_mapping, clone_pal = clone_pal, show_clone_label = show_clone_label,
      idx = sample_label_idx, show_legend = show_legend, show_library_label = show_library_label,
      str_to_remove = str_to_remove
    ),
    heatmap_legend_param = leg_params,
    top_annotation = make_top_annotation_gain(copynumber,
      cutoff = cutoff, maxf = maxf,
      plotfrequency = plotfrequency, plotcol = plotcol, SV = SV
    ),
    use_raster = TRUE,
    raster_quality = 10,
    ...
  )
  return(copynumber_hm)
}

getSVlegend <- function(include = NULL) {
  svs <- SV_colors[include]
  SV <- list(
    ComplexHeatmap::Legend(
      labels = names(svs), title = "Rearrangement type", type = "points", pch = 16,
      legend_gp = grid::gpar(col = as.vector(svs))
    )
  )
  return(SV)
}

#' Heatmap plot
#'
#' Plot a heatmap where rows are cells, columns are genome coordinates and colours map to (allele-specific) copy-number states
#'
#' @param cn Either a hscn object or a single cell allele specific copy number dataframe with the following columns: `cell_id`, `chr`, `start`, `end`, `state`, `copy`
#' @param tree Tree in newick format to plot alongside the heatmap, default = NULL
#' @param cluster data.frame assigning cells to clusters, needs the following columns `cell_id`, `clone_id` default = NULL
#' @param normalize_ploidy Normalize ploidy of all cells to 2
#' @param normalize_tree default = FALSE
#' @param branch_length scales branch lengths to this size, default = 2
#' @param spacer_cols number of empty columns between chromosomes, default = 20
#' @param plottree Binary value of whether to plot tree or not, default = TRUE
#' @param plotcol Which column to colour the heatmap by, should be one of "state", "state_BAF", "state_phase", "state_AS", "state_min", "copy", "BAF", "Min", "Maj"
#' @param reordercluster Reorder the cells according to cluster if no tree is specified
#' @param pctcells Minimum size of cluster in terms of perecentage of cells in umap clustering
#' @param library_mapping Named vector mapping library names to labels for legend
#' @param clone_pal pallette to colour clusters by
#' @param sample_label_idx default = 1
#' @param fillna Smooth over NA values, default = TRUE
#' @param frequencycutoff default = 2
#' @param maxf Max frequency when plotting the frequency track, default = NULL infers this from the data
#' @param plotfrequency Plot the frequency track of gains and losses across the genome
#' @param show_legend plot legend or not, boolean
#' @param show_library_label show library label or not, boolean
#' @param show_clone_label show clone label or not, boolean
#' @param umapmetric metric to use in umap dimensionality reduction if no clusters are specified
#' @param chrlabels include chromosome labels or not, boolean
#' @param SV sv data frame
#' @param seed seed for UMAP
#' @param nticks number of ticks in x-axis label when plotting a single chromosome
#' @param fillgenome fill in any missing bins and add NA to centromeric regions
#' @param na_col colour of NA values
#' @param linkheight height of x-axis ticks
#' @param newlegendname overwrite default legend name
#' @param str_to_remove string to remove from cell_id's when plotting labels
#'
#' If clusters are set to NULL then the function will compute clusters using UMAP and HDBSCAN.
#'
#' #' \dontrun{
#' data("haplotypes")
#' data("CNbins")
#' haplotypes <- format_haplotypes_dlp(haplotypes, CNbins)
#' hscn <- callHaplotypeSpecificCN(CNbins, haplotypes, likelihood = "binomial")
#' plotHeatmap(hscn)
#' }
#'
#' @export
plotHeatmap <- function(cn,
                        tree = NULL,
                        clusters = NULL,
                        normalize_ploidy = FALSE,
                        normalize_tree = FALSE,
                        branch_length = 1,
                        spacer_cols = 20,
                        plottree = TRUE,
                        plotcol = "state",
                        reorderclusters = FALSE,
                        pctcells = 0.05,
                        library_mapping = NULL,
                        clone_pal = NULL,
                        sample_label_idx = 1,
                        fillna = TRUE,
                        frequencycutoff = 2,
                        maxf = NULL,
                        plotfrequency = FALSE,
                        show_legend = TRUE,
                        show_library_label = TRUE,
                        show_clone_label = TRUE,
                        widenarm = FALSE,
                        umapmetric = "euclidean",
                        chrlabels = TRUE,
                        SV = NULL,
                        seed = NULL,
                        nticks = 4,
                        fillgenome = FALSE,
                        annotation_height = NULL, 
                        annofontsize = 14,
                        na_col = "white",
                        linkheight = 5,
                        newlegendname = NULL,
                        str_to_remove = NULL,
                        ...) {
  if (is.hscn(cn) | is.ascn(cn)) {
    CNbins <- cn$data
  } else {
    CNbins <- cn
  }

  if (widenarm == TRUE) {
    dlpbinsarm <- dlpbins %>%
      dplyr::mutate(arm = coord_to_arm(chr, start), chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(arm = ifelse(chrarm %in% unique(CNbins$chrarm), arm, "")) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      as.data.table()

    dlpbinsarm <- data.table::rbindlist(lapply(
      unique(CNbins$cell_id),
      function(i) {
        cbind(dlpbinsarm,
          cell_id = i
        )
      }
    )) %>%
      data.table::setkey("chr", "arm", "chrarm", "start", "end")

    CNbinst <- setkey(as.data.table(CNbins %>% dplyr::select(-start, -end)), "chr", "arm", "chrarm")
    CNbins <- dlpbinsarm[CNbinst, on = c("chr", "chrarm", "arm", "cell_id")] %>%
      .[!is.na(cell_id)] %>%
      orderdf(.)
  }

  if (!plotcol %in% c("state", "state_BAF", "state_phase", "state_AS", "state_min", "copy", "BAF", "Min", "Maj")) {
    stop(paste0("Column name - ", plotcol, " not available for plotting, please use one of state, copy, BAF, state_BAF, state_phase, state_AS, Min or Maj"))
  }

  if (!plotcol %in% names(CNbins)) {
    stop(paste0("Column name - ", plotcol, " not in CNbins data frame..."))
  }

  if (plotcol == "state") {
    colvals <- cn_colours
    legendname <- "Copy Number"
  }

  if (plotcol == "Min") {
    colvals <- cn_colours
    legendname <- "Copy Number\nAllele B"
  }

  if (plotcol == "Maj") {
    colvals <- cn_colours
    legendname <- "Copy Number\nAllele A"
  }

  if (plotcol == "state_BAF") {
    colvals <- cn_colours_bafstate
    legendname <- "Allelic Imbalance"
  }

  if (plotcol == "BAF") {
    colvals <- circlize::colorRamp2(c(0, 0.5, 1), c(scCNphase_colors["A-Hom"], scCNphase_colors["Balanced"], scCNphase_colors["B-Hom"]))
    legendname <- "Allelic Imbalance"
  }

  if (plotcol == "copy") {
    colvals <- circlize::colorRamp2(seq(0, 11, 1), scCN_colors)
    legendname <- "Copy"
  }

  if (plotcol == "state_AS") {
    colvals <- cn_colours_loh
    legendname <- "Allele Specific Copy Number"
  }

  if (plotcol == "state_min") {
    colvals <- cn_colours_minorallele
    legendname <- "Minor Allele Copy Number"
  }

  if (plotcol == "state_phase") {
    colvals <- cn_colours_phase
    legendname <- "Allelic Imbalance"
  }
  
  if (!is.null(newlegendname)){
    legendname <- newlegendname
  }

  ncells <- length(unique(CNbins$cell_id))

  if (is.null(clusters)) {
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
  } else {
    ordered_cell_ids <- paste0(clusters$cell_id)
  }

  if (is.null(tree) & is.null(clusters)) {
    message("No tree or cluster information provided, clustering using HDBSCAN")
    clustering_results <- umap_clustering(CNbins,
      minPts = max(round(pctcells * ncells), 2),
      field = "copy",
      umapmetric = umapmetric,
      seed = seed
    )
    tree <- clustering_results$tree
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)

    clusters <- clustering_results$clustering %>%
      dplyr::select(cell_id, clone_id)
  }

  if (!is.null(clusters)) {
    cells_clusters <- length(unique(clusters$cell_id))
    cells_data <- length(unique(CNbins$cell_id))
    if (cells_data != cells_clusters){
      warning("Number of cells in clusters dataframe !=  number of cells in the bins data!")
    }
    if (!"clone_id" %in% names(clusters)) {
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    if (reorderclusters == TRUE) {
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  if (plottree == TRUE) {
    if (normalize_tree == T) {
      tree <- format_tree(tree, branch_length)
    }

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data

    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
  }

  if (!is.null(clusters)) {
    if (!"clone_id" %in% names(clusters)) {
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    else if (reorderclusters == TRUE & !is.null(tree)) {
      message("Reorder clusters dataframe according to clones using tree")
      if (normalize_tree == T) {
        tree <- format_tree(tree, branch_length)
      }

      tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
      tree_plot_dat <- tree_ggplot$data

      message("Creating tree...")
      tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
      ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
    }
    else if (reorderclusters == TRUE & is.null(tree)) {
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  message("Creating copy number heatmap...")
  if (fillgenome) {
    copynumber <- createCNmatrix(CNbins, field = plotcol, wholegenome = TRUE,
                                 fillnaplot = fillna, centromere = FALSE)
  } else {
    copynumber <- createCNmatrix(CNbins, field = plotcol, fillnaplot = fillna)
  }
  if (normalize_ploidy == T) {
    message("Normalizing ploidy for each cell to 2")
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  copynumber <- format_copynumber(copynumber,
    ordered_cell_ids,
    spacer_cols = spacer_cols,
    plotcol = plotcol
  )
  clones_formatted <- format_clones(as.data.frame(clusters), ordered_cell_ids)
  if (!is.null(clone_pal)) {
    clones_idx <- dplyr::distinct(clones_formatted, clone_id, clone_label)
    clone_pal <- clone_pal[clones_idx$clone_id]
    names(clone_pal) <- clones_idx$clone_label
  }
  copynumber_hm <- make_copynumber_heatmap(copynumber,
    clones_formatted,
    colvals = colvals,
    legendname = legendname,
    library_mapping = library_mapping,
    clone_pal = clone_pal,
    sample_label_idx = sample_label_idx,
    cutoff = frequencycutoff,
    maxf = maxf,
    plotcol = plotcol,
    plotfrequency = plotfrequency,
    show_legend = show_legend,
    show_library_label = show_library_label,
    show_clone_label = show_clone_label,
    chrlabels = chrlabels,
    SV = SV,
    nticks = nticks,
    annotation_height = annotation_height, 
    annofontsize = annofontsize,
    na_col = na_col,
    linkheight = linkheight,
    str_to_remove = str_to_remove,
    ...
  )
  if (plottree == TRUE) {
    h <- tree_hm + copynumber_hm
  } else {
    h <- copynumber_hm
  }

  return(h)
}

#' @export
createSNVmatrix <- function(SNVs, allcells = NULL, field = "VAF") {
  if ("clone_id" %in% names(SNVs)) {
    SNVs$cell_id <- SNVs$clone_id
  }

  snvmatrix <- SNVs %>%
    .[, mutid := paste(chr, as.integer(start), ref, alt, sep = "_")] %>%
    data.table::dcast(., cell_id ~ mutid, value.var = field, fill = 0) %>%
    as.data.frame()

  if (!is.null(allcells)) {
    message("Adding blank rows to cells that have no mutations...")
    missingcells <- data.frame(cell_id = clones$cell_id[!clones$cell_id %in% dfmuts$cell_id])
    dfmuts <- dplyr::bind_rows(dfmuts, missingcells)
  }


  rownames(snvmatrix) <- snvmatrix$cell_id
  snvmatrix <- subset(snvmatrix, select = -cell_id)

  # sort by number of clones with mutation
  snvmatrix <- snvmatrix[, names(sort(colSums(!is.na(snvmatrix)), decreasing = T))]

  return(snvmatrix)
}

#' @export
plotSNVHeatmap <- function(SNVs,
                           tree,
                           clusters = NULL,
                           field = "VAF",
                           reorderclusters = FALSE,
                           clone_pal = NULL,
                           show_legend = TRUE,
                           library_mapping = NULL,
                           show_library_label = TRUE,
                           show_clone_label = TRUE,
                           plottree = TRUE,
                           mymaxcol = "firebrick4",
                           sample_label_idx = 1,
                           nsample = 10000,
                           clustercols = FALSE,
                           str_to_remove = NULL) {
  muts <- createSNVmatrix(SNVs)

  if (is.null(clusters)) {
    clusters <- data.frame(cell_id = tree$tip.label, clone_id = "0")
  }

  tree_ggplot <- make_tree_ggplot(tree, clusters, clone_pal = clone_pal)
  tree_plot_dat <- tree_ggplot$data

  message("Creating tree...")
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
  ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)

  muts <- muts[ordered_cell_ids, ]

  cols <- circlize::colorRamp2(c(0, 1), c("white", mymaxcol))

  clones_formatted <- format_clones(as.data.frame(clusters), ordered_cell_ids)

  muts <- as.matrix(muts)

  if (dim(muts)[2] > 10000) {
    message(paste0("Sampling ", nsample, " mutations..."))
    muts <- muts[, sample(ncol(muts), size = nsample), drop = FALSE]
  }

  snv_hm <- ComplexHeatmap::Heatmap(
    name = field,
    muts,
    col = cols,
    na_col = "grey20",
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = clustercols,
    show_column_names = FALSE,
    # bottom_annotation=make_bottom_annot(copynumber),
    left_annotation = make_left_annot(muts, clones_formatted,
      library_mapping = library_mapping, clone_pal = clone_pal, show_clone_label = show_clone_label,
      idx = sample_label_idx, show_legend = show_legend, show_library_label = show_library_label,
      str_to_remove = strstr_to_remove
    ),
    use_raster = TRUE,
    # top_annotation = HeatmapAnnotation(df = mutgroups,col = list(MutationGroup = colpal)),
    heatmap_legend_param = list(nrow = 4)
  )

  if (plottree == TRUE) {
    h <- tree_hm + snv_hm
  } else {
    h <- snv_hm
  }

  return(h)
}

plotHeatmapQC <- function(cn,
                          tree = NULL,
                          clusters = NULL,
                          normalize_ploidy = FALSE,
                          normalize_tree = FALSE,
                          branch_length = 1,
                          spacer_cols = 20,
                          plottree = TRUE,
                          plotcol = "state",
                          reorderclusters = FALSE,
                          pctcells = 0.05,
                          library_mapping = NULL,
                          clone_pal = NULL,
                          sample_label_idx = 1,
                          plotfrequency = FALSE,
                          show_legend = TRUE,
                          show_library_label = TRUE,
                          ...) {
  CNbins <- cn$data
  p <- cn$phasing

  if (is.null(clusters)) {
    clusters <- p$cl$clustering
    tree <- p$cl$tree
  }

  if ("chrarm" %in% names(p$prop)) {
    arm <- TRUE
  } else {
    arm <- FALSE
  }

  if (arm == FALSE) {
    CNbins <- dplyr::left_join(CNbins, p$cl$clustering, by = "cell_id")
    CNbins <- dplyr::left_join(CNbins, p$prop)
  } else {
    CNbins$chrarm <- paste0(CNbins$chr, schnapps:::coord_to_arm(CNbins$chr, CNbins$start))
    CNbins <- dplyr::left_join(CNbins, p$cl$clustering, by = "cell_id")
    CNbins <- dplyr::left_join(CNbins, p$prop)
  }

  CNbins <- CNbins %>%
    dplyr::mutate(state_BAF = ifelse(is.na(propA), "-1", state_BAF))

  plotcol <- "state_BAF"
  colvals <- schnapps:::cn_colours_bafstate
  colvals[["-1"]] <- "gray90"
  legendname <- "Allelic Imbalance"

  ncells <- length(unique(CNbins$cell_id))

  if (is.null(clusters)) {
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
  } else {
    ordered_cell_ids <- paste0(clusters$cell_id)
  }

  if (is.null(tree) & is.null(clusters)) {
    message("No tree or cluster information provided, clustering using HDBSCAN")
    clustering_results <- umap_clustering(CNbins, minPts = max(round(pctcells * ncells), 2), field = "copy")
    tree <- clustering_results$tree
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)

    clusters <- clustering_results$clustering %>%
      dplyr::select(cell_id, clone_id)
  }

  if (!is.null(clusters)) {
    if (!"clone_id" %in% names(clusters)) {
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    if (reorderclusters == TRUE) {
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  if (plottree == TRUE) {
    if (normalize_tree == T) {
      tree <- format_tree(tree, branch_length)
    }

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data

    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
  }

  message("Creating copy number heatmap...")
  copynumber <- createCNmatrix(CNbins, field = plotcol, fillna = TRUE)
  if (normalize_ploidy == T) {
    message("Normalizing ploidy for each cell to 2")
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  copynumber <- format_copynumber(copynumber, ordered_cell_ids, spacer_cols = spacer_cols)
  clones_formatted <- format_clones(as.data.frame(clusters), ordered_cell_ids)
  if (!is.null(clone_pal)) {
    clones_idx <- dplyr::distinct(clones_formatted, clone_id, clone_label)
    clone_pal <- clone_pal[clones_idx$clone_id]
    names(clone_pal) <- clones_idx$clone_label
  }
  copynumber_hm <- make_copynumber_heatmap(copynumber,
    clones_formatted,
    colvals = colvals,
    legendname = legendname,
    library_mapping = library_mapping,
    clone_pal = clone_pal,
    plotfrequency = plotfrequency,
    show_legend = show_legend,
    show_library_label = show_library_label,
    sample_label_idx = sample_label_idx
  )
  if (plottree == TRUE) {
    h <- tree_hm + copynumber_hm
  } else {
    h <- copynumber_hm
  }

  return(h)
}
