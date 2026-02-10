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

#' @export
make_copynumber_legend <- function(font_size = 12, ncolcn = 2, ncolas = 1, gainloss = FALSE, cnonly = FALSE, cntitle = "Copy\nNumber", hscntitle = "HSCN State", gainlosstitle = "D", ...) {
  cn_lgd <- ComplexHeatmap::Legend(
    title = cntitle,
    labels = stringr::str_remove(names(scCN_colors), "CN"),
    legend_gp = grid::gpar(fill = as.vector(scCN_colors)),
    labels_gp = grid::gpar(fontsize = font_size),
    title_gp = grid::gpar(fontsize = font_size),
    title_gap = grid::unit(1, "mm"),
    grid_height = grid::unit(3, "mm"), grid_width = grid::unit(2.5, "mm"),
    ncol = ncolcn
  )
  
  hscn_lgd <- ComplexHeatmap::Legend(
    title = hscntitle,
    labels = names(scCNphase_colors),
    legend_gp = grid::gpar(fill = as.vector(scCNphase_colors)),
    labels_gp = grid::gpar(fontsize = font_size),
    title_gp = grid::gpar(fontsize = font_size),
    title_gap = grid::unit(1, "mm"),
    grid_height = grid::unit(3, "mm"), grid_width = grid::unit(2.5, "mm"),
    ncol = ncolas
  )
  
  gain_loss_lgd <- ComplexHeatmap::Legend(
    title = gainlosstitle,
    labels = c("Gain", "Loss"),
    legend_gp = grid::gpar(fill = c("#E34A33", "#3182BD")),
    labels_gp = grid::gpar(fontsize = font_size),
    title_gp = grid::gpar(fontsize = font_size),
    title_gap = grid::unit(1, "mm"),
    grid_height = grid::unit(3, "mm"), grid_width = grid::unit(2.5, "mm"),
    ncol = 1
  )
  
  if (gainloss){
    lgd <- ComplexHeatmap::packLegend(
      cn_lgd, hscn_lgd, gain_loss_lgd,
      row_gap = grid::unit(4, "mm"),
      column_gap = grid::unit(4, "mm"),
      ...
    )
  } else {
    lgd <- ComplexHeatmap::packLegend(
      cn_lgd, hscn_lgd,
      row_gap = grid::unit(4, "mm"),
      column_gap = grid::unit(4, "mm"),
      ...
    )
  }
  
  if (cnonly){
    lgd <- cn_lgd
  }
  
  grid::grid.grabExpr(ComplexHeatmap::draw(lgd))
}

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

make_tree_ggplot <- function(tree, clones, clone_pal = NULL, ladderize = TRUE) {
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

  p <- ggtree::ggtree(tree, tree_aes, size = 0.25, ladderize = ladderize) +
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

  if (plotcol %in% c("BAF", "copy", "other")) {
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

make_corrupt_tree_heatmap <- function(tree_ggplot, tree_width, ...) {
  tree_annot_func <- ComplexHeatmap::AnnotationFunction(
    fun = function(index) {
      grid::pushViewport(grid::viewport(height = 1))

      # Avoid hard-coded grob indices: grob order is not stable across ggplot2/ggtree versions.
      tree_grob <- ggplot2::ggplotGrob(tree_ggplot)
      panel_idx <- which(tree_grob$layout$name == "panel")
      if (length(panel_idx) == 0) {
        panel_idx <- which(grepl("^panel", tree_grob$layout$name))
      }

      if (length(panel_idx) >= 1) {
        grid::grid.draw(tree_grob$grobs[[panel_idx[[1]]]])
      } else {
        # Fallback: draw the full grob if we can't reliably locate the panel.
        grid::grid.draw(tree_grob)
      }

      grid::popViewport()
    },
    var_import = list(tree_ggplot = tree_ggplot),
    width = grid::unit(tree_width, "cm"),
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

#' @export
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

make_left_annot_generic <- function(dfanno,
                                   palettes = NULL,
                                   show_legend = TRUE,
                                   annofontsize = 14,
                                   anno_width = 0.4) {
  # Convert to data.frame if needed
  if (inherits(dfanno, c("data.table", "tbl_df", "tbl"))) {
    message("Converting dfanno from ", class(dfanno)[1], " to data.frame")
    dfanno <- as.data.frame(dfanno)
  }
  if (!is.data.frame(dfanno)) {
    warning("dfanno was converted to data.frame from ", class(dfanno)[1])
    dfanno <- as.data.frame(dfanno)
  }
  
  # Check if cell_id column exists
  if (!"cell_id" %in% colnames(dfanno)) {
    stop("dfanno must contain a 'cell_id' column")
  }
  
  # Check for NA values in dfanno
  na_cols <- colnames(dfanno)[sapply(dfanno, function(x) any(is.na(x)))]
  if (length(na_cols) > 0) {
    stop(paste0("dfanno contains NA values in the following column(s): ", 
                paste(na_cols, collapse = ", "), 
                ". Please remove or impute NA values before calling this function."))
  }
  
  # Get annotation columns (all except cell_id)
  anno_cols <- setdiff(colnames(dfanno), "cell_id")
  if (length(anno_cols) == 0) {
    stop("dfanno must contain at least one annotation column besides cell_id")
  }
  
  # Detect continuous vs discrete columns
  continuous_cols <- character(0)
  discrete_cols <- character(0)
  col_types <- list()
  
  for (col in anno_cols) {
    col_data <- dfanno[[col]]
    is_numeric <- is.numeric(col_data) || is.integer(col_data)
    n_unique <- length(unique(col_data))
    
    if (is_numeric && n_unique > 10) {
      continuous_cols <- c(continuous_cols, col)
      col_types[[col]] <- "continuous"
    } else {
      discrete_cols <- c(discrete_cols, col)
      col_types[[col]] <- "discrete"
    }
  }
  
  # Default palette rotation if none specified
  if (is.null(palettes)) {
    default_palettes <- c("Set2", "Set1", "Set3", "Paired", "Dark2", "Accent")
  } else {
    default_palettes <- palettes
  }
  
  # Create annotation colors for discrete columns
  annot_colours <- list()
  for (i in seq_along(discrete_cols)) {
    col <- discrete_cols[i]
    # Cycle through palettes
    palette_idx <- ((i-1) %% length(default_palettes)) + 1
    current_palette <- default_palettes[palette_idx]
    
    # Get unique values for this column
    levels <- gtools::mixedsort(unique(dfanno[[col]]))
    # Create palette for this annotation
    annot_colours[[col]] <- make_discrete_palette(current_palette, levels)
  }
  
  # Create color mappings for continuous columns
  continuous_col_funs <- list()
  if (length(continuous_cols) > 0) {
    # Check if viridis is available, otherwise use a fallback
    if (requireNamespace("viridis", quietly = TRUE)) {
      viridis_colors <- viridis::viridis(100)
    } else {
      # Fallback: use viridis-like gradient via colorRampPalette
      viridis_colors <- grDevices::colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE725"))(100)
    }
    
    for (col in continuous_cols) {
      col_data <- dfanno[[col]]
      col_range <- range(col_data, na.rm = TRUE)
      
      # Create color mapping function and store it
      continuous_col_funs[[col]] <- circlize::colorRamp2(
        seq(col_range[1], col_range[2], length.out = 100),
        viridis_colors
      )
      
      # Add to annot_colours for use with df parameter
      annot_colours[[col]] <- continuous_col_funs[[col]]
    }
  }
  
  # Build annotation using df parameter (auto-generates legends for all columns)
  # This works for both discrete (named color vectors) and continuous (color functions)
  legend_params <- list()
  for (col in anno_cols) {
    if (col_types[[col]] == "discrete") {
      legend_params[[col]] <- list(
        nrow = 3,
        direction = "horizontal",
        labels_gp = grid::gpar(fontsize = annofontsize-1),
        title_gp = grid::gpar(fontsize = annofontsize-1),
        legend_gp = grid::gpar(fontsize = annofontsize-1)
      )
    } else {
      # Continuous columns get simpler legend params
      legend_params[[col]] <- list(
        labels_gp = grid::gpar(fontsize = annofontsize-1),
        title_gp = grid::gpar(fontsize = annofontsize-1),
        legend_gp = grid::gpar(fontsize = annofontsize-1)
      )
    }
  }
  
  # Create the annotation object using df parameter
  # This approach auto-generates legends for both discrete and continuous columns
  left_annot <- ComplexHeatmap::HeatmapAnnotation(
    df = as.data.frame(dfanno[, anno_cols, drop = FALSE]),
    col = annot_colours,
    which = "row",
    show_legend = show_legend,
    annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
    simple_anno_size = grid::unit(anno_width, "cm"),
    annotation_legend_param = legend_params
  )
  
  return(left_annot)
}

make_left_annot <- function(copynumber,
                            clones,
                            library_mapping = NULL,
                            show_library_label = TRUE,
                            show_clone_label = TRUE,
                            show_clone_text = TRUE,
                            clone_pal = NULL,
                            idx = 1,
                            show_legend = TRUE,
                            annofontsize = 14,
                            anno_width = 0.4,
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
  library_levels <- gtools::mixedsort(unlist(unique(library_labels)))
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
        gp = grid::gpar(fontsize = annofontsize - 1),
        just = c("centre", "centre")
      )
    }

    clone_legend_rows <- 3
    if (length(clone_levels) > 3) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 2))
    }

    if (show_library_label == TRUE & show_clone_label == TRUE & show_clone_text == TRUE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Cluster = clones$clone_label, 
        clone_label = clone_label_generator,
        Sample = library_labels,
        col = annot_colours, show_annotation_name = c(TRUE, FALSE, TRUE),
        which = "row", annotation_width = grid::unit(rep(anno_width, 3), "cm"),
        annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
        annotation_legend_param = list(
          Cluster = list(nrow = clone_legend_rows, direction = "horizontal"),
          Sample = list(nrow = library_legend_rows, direction = "horizontal"),
          labels_gp = grid::gpar(fontsize = annofontsize-1),
          legend_gp = grid::gpar(fontsize = annofontsize-1),
          title_gp = grid::gpar(fontsize = annofontsize-1)
        ),
        show_legend = show_legend
      )
    } else if (show_library_label == FALSE & show_clone_label == TRUE & show_clone_text == TRUE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Cluster = clones$clone_label, 
        clone_label = clone_label_generator,
        col = annot_colours, show_annotation_name = c(TRUE, FALSE),
        which = "row", annotation_width = grid::unit(rep(anno_width, 2), "cm"),
        annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
        annotation_legend_param = list(
          Cluster = list(nrow = clone_legend_rows),
          labels_gp = grid::gpar(fontsize = annofontsize-1),
          legend_gp = grid::gpar(fontsize = annofontsize-1),
          title_gp = grid::gpar(fontsize = annofontsize-1)
        ),
        show_legend = show_legend
      ) 
      } else if (show_library_label == TRUE & show_clone_label == TRUE & show_clone_text == FALSE){
        left_annot <- ComplexHeatmap::HeatmapAnnotation(
          Cluster = clones$clone_label, 
          Sample = library_labels,
          col = annot_colours, show_annotation_name = c(TRUE, FALSE, TRUE),
          which = "row", annotation_width = grid::unit(rep(anno_width, 3), "cm"),
          annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
          annotation_legend_param = list(
            Cluster = list(nrow = clone_legend_rows, direction = "horizontal"),
            Sample = list(nrow = library_legend_rows, direction = "horizontal"),
            labels_gp = grid::gpar(fontsize = annofontsize-1),
            legend_gp = grid::gpar(fontsize = annofontsize-1),
            title_gp = grid::gpar(fontsize = annofontsize-1)
          ),
          show_legend = show_legend
        )
      } else if (show_library_label == TRUE & show_clone_label == FALSE) {
        left_annot <- ComplexHeatmap::HeatmapAnnotation(
          Sample = library_labels, col = annot_colours,
          which = "row", simple_anno_size = grid::unit(anno_width, "cm"),
          annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
          annotation_legend_param = list(
            Sample = list(nrow = library_legend_rows),
            labels_gp = grid::gpar(fontsize = annofontsize-1),
            legend_gp = grid::gpar(fontsize = annofontsize-1),
            title_gp = grid::gpar(fontsize = annofontsize-1)
          ),
          show_legend = show_legend
        )
      } else if (show_library_label == FALSE & show_clone_label == TRUE & show_clone_text == FALSE) {
        left_annot <- ComplexHeatmap::HeatmapAnnotation(
          Cluster = clones$clone_label, 
          col = annot_colours, show_annotation_name = c(TRUE, FALSE),
          which = "row", annotation_width = grid::unit(rep(anno_width, 2), "cm"),
          simple_anno_size = grid::unit(anno_width * 2, "cm"),
          annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
          annotation_legend_param = list(
            Cluster = list(nrow = clone_legend_rows),
            labels_gp = grid::gpar(fontsize = annofontsize-1),
            legend_gp = grid::gpar(fontsize = annofontsize-1),
            title_gp = grid::gpar(fontsize = annofontsize-1)
          ),
          show_legend = show_legend
        ) 
      }
  } else {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Sample = library_labels, col = annot_colours,
      which = "row", simple_anno_size = grid::unit(0.4, "cm"),
      annotation_name_gp = grid::gpar(fontsize = annofontsize - 1),
      annotation_legend_param = list(
        Cluster = list(nrow = clone_legend_rows),
        labels_gp = grid::gpar(fontsize = annofontsize-1),
        legend_gp = grid::gpar(fontsize = annofontsize-1),
        title_gp = grid::gpar(fontsize = annofontsize-1)
      ),
      show_legend = show_legend
    )
  }

  return(left_annot)
}

make_top_annotsnv <- function(mutgroups) {

}

get_genomecoords_label_pos <- function(copynumber, Mb = TRUE, nticks = 3) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  binwidth <- as.numeric(strsplit(colnames(copynumber)[1], ":")[[1]][3]) -
    as.numeric(strsplit(colnames(copynumber)[1], ":")[[1]][2]) + 1
  chromfreq <- table(chroms)
  uniq_chroms <- names(chromfreq)[chromfreq > 1]
  uniq_chroms <- uniq_chroms[stringr::str_detect(uniq_chroms, "V", negate = TRUE)]
  chromvec <- lapply(colnames(copynumber), function(x) strsplit(x, ":")[[1]][1])
  binvec <- lapply(colnames(copynumber), function(x) (as.numeric(strsplit(x, ":")[[1]][2]) - 1) / 1e6)
  binvec <- split(unlist(binvec), unlist(chromvec))
  for (chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    binvec_ <- binvec[[chrom]]
    max_idx <- round(0.9 * max(binvec_))
    nwidth <- 10 * round((max(binvec_) / nticks) / 10)
    mypos <- nwidth
    for (ticks in 1:nticks){
      if (mypos > max_idx){
        next
      }
      if (Mb){
        lab <- paste0(mypos, "Mb")
      } else {
        lab <- paste0(mypos)
      }
      chrom_label_pos[[lab]] <- which.min(abs(binvec_ - mypos))
      mypos <- mypos + nwidth
    }
  }
  return(chrom_label_pos)
}

get_chrom_label_pos <- function(copynumber, Mb = TRUE, nticks = 3) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  chromfreq <- table(chroms)
  uniq_chroms <- names(chromfreq)[chromfreq > 1]
  uniq_chroms <- uniq_chroms[stringr::str_detect(uniq_chroms, "V", negate = TRUE)]

  if (length(uniq_chroms) == 1){
   chrom_label_pos <- get_genomecoords_label_pos(copynumber, Mb = Mb, nticks = nticks)
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
                              Mb = TRUE,
                              labeladjust = -1,
                              annotation_height = NULL, 
                              annofontsize = 14,
                              linkheight = 1) {
  if (chrlabels[1] == FALSE) {
    return(NULL)
  } else if (chrlabels[1] == TRUE) {
    chrom_label_pos <- get_chrom_label_pos(copynumber, Mb = Mb, nticks = nticks)
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels = ComplexHeatmap::anno_mark(
      at = as.vector(unlist(chrom_label_pos)),
      labels = names(chrom_label_pos),
      link_height = grid::unit(linkheight, "mm"),
      labels_gp = grid::gpar(fontsize = annofontsize),
      side = "bottom",
      padding = grid::unit(labeladjust, "mm"), 
      extend = 0.01, 
      labels_rot = 0
    ), show_annotation_name = FALSE,
    annotation_height = annotation_height)
  } else {
    chrom_label_pos <- get_chrom_label_pos(copynumber)
    chrom_label_pos <- chrom_label_pos[chrlabels]
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels = ComplexHeatmap::anno_mark(
      at = as.vector(unlist(chrom_label_pos)),
      labels = names(chrom_label_pos),
      link_height = grid::unit(linkheight, "mm"),
      side = "bottom",
      labels_gp = grid::gpar(fontsize = annofontsize),
      padding = grid::unit(labeladjust, "mm"), extend = 0.01, labels_rot = 0
    ), show_annotation_name = FALSE,
    annotation_height = annotation_height)
  }
  return(bottom_annot)
}

#' Create ideogram annotation for heatmap
#'
#' Creates a ComplexHeatmap annotation showing chromosome ideogram (cytobands)
#' at the bottom of the heatmap.
#'
#' @param copynumber Copy number matrix with column names in format "chr:start:end"
#' @param genome Genome assembly ("hg19" or "hg38"). Default is "hg19".
#' @param ideogram_height Height of the ideogram annotation in cm. Default is 0.3.
#'
#' @return A ComplexHeatmap HeatmapAnnotation object
#'
#' @keywords internal
make_ideogram_annotation <- function(copynumber,
                                     genome = "hg19",
                                     ideogram_height = 0.3) {
  # Load cytoband data
  data("cytoband_map", envir = environment())
  cytobands <- cytoband_map[[genome]]

  if (is.null(cytobands) || nrow(cytobands) == 0) {
    warning("No cytoband data found for genome: ", genome)
    return(NULL)
  }

  # Parse column names to get coordinates (format: "chr:start:end" or spacer columns "V*")
  col_names <- colnames(copynumber)
  n_cols <- length(col_names)

  # Initialize stain type vector for all columns (use "spacer" for gaps)
  stain_types <- rep("spacer", n_cols)

  # Process each column
  for (i in seq_len(n_cols)) {
    col_name <- col_names[i]

    # Skip spacer columns (contain "V" pattern from space_copynumber_columns)
    if (grepl("^V[0-9]+$", col_name) || is.na(col_name)) {
      next
    }

    # Parse coordinate from column name
    parts <- strsplit(col_name, ":")[[1]]
    if (length(parts) != 3) {
      next
    }

    chr <- parts[1]
    bin_start <- as.numeric(parts[2])
    bin_end <- as.numeric(parts[3])
    bin_mid <- (bin_start + bin_end) / 2

    # Find matching cytoband (use midpoint of bin)
    chr_with_prefix <- paste0("chr", chr)
    matching_band <- cytobands[V1 == chr_with_prefix & V2 <= bin_mid & V3 > bin_mid]

    if (nrow(matching_band) > 0) {
      stain <- matching_band$V5[1]
      if (stain %in% names(cyto_colors)) {
        stain_types[i] <- stain
      }
    }
  }

  # Create color mapping with spacer as white
  all_colors <- c(cyto_colors, spacer = "white")

  # Create the annotation using a custom annotation function for colored bars
  ideogram_annot <- ComplexHeatmap::HeatmapAnnotation(
    ideogram = ComplexHeatmap::anno_simple(
      stain_types,
      col = all_colors,
      height = grid::unit(ideogram_height, "cm")
    ),
    which = "column",
    show_annotation_name = FALSE,
    show_legend = FALSE
  )

  return(ideogram_annot)
}

#' Find column positions for gene annotations
#'
#' Maps gene names to their corresponding column indices in the copynumber matrix.
#' Genes are matched to bins where the gene start position falls within the bin boundaries.
#'
#' @param copynumber Copy number matrix with column names in format "chr:start:end"
#' @param gene_annotations Character vector of gene names to annotate
#' @param genome Genome assembly ("hg19" or "hg38"). Default is "hg19".
#' @param gene_label_sep Separator when multiple genes fall in the same bin. Default is "/".
#'
#' @return A list with components:
#'   \item{at}{Integer vector of column indices}
#'   \item{labels}{Character vector of gene labels}
#'   Returns NULL if no genes are found in the matrix columns.
#'
#' @keywords internal
find_gene_bin_positions <- function(copynumber,
                                    gene_annotations,
                                    genome = "hg19",
                                    gene_label_sep = "/") {
  # Load gene locations data

  data("gene_locations", envir = environment())
  gene_data <- gene_locations[[genome]]

  if (is.null(gene_data) || nrow(gene_data) == 0) {
    warning("No gene location data found for genome: ", genome)
    return(NULL)
  }

  # Parse column names to get bin coordinates
  col_names <- colnames(copynumber)
  n_cols <- length(col_names)

  # Build a data frame of non-spacer columns with their indices
  col_info <- data.frame(
    col_idx = integer(0),
    chr = character(0),
    start = numeric(0),
    end = numeric(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_cols)) {
    col_name <- col_names[i]

    # Skip spacer columns (V1, V2, etc.)
    if (grepl("^V[0-9]+$", col_name) || is.na(col_name)) {
      next
    }

    # Parse coordinate from column name (format: "chr:start:end")
    parts <- strsplit(col_name, ":")[[1]]
    if (length(parts) != 3) {
      next
    }

    col_info <- rbind(col_info, data.frame(
      col_idx = i,
      chr = parts[1],
      start = as.numeric(parts[2]),
      end = as.numeric(parts[3]),
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(col_info) == 0) {
    warning("No valid genomic bins found in copynumber matrix")
    return(NULL)
  }

  # Find bin positions for each gene
  gene_positions <- list()  # col_idx -> list of gene names
  genes_not_found <- character(0)
  genes_not_in_bins <- character(0)

  for (gene in gene_annotations) {
    # Find gene in gene_data
    gene_row <- gene_data[gene_data$ensembl_gene_symbol == gene, ]

    if (nrow(gene_row) == 0) {
      genes_not_found <- c(genes_not_found, gene)
      next
    }

    # Use first occurrence if multiple (some genes have multiple loci)
    gene_chr <- gene_row$chr[1]
    gene_start <- gene_row$start[1]

    # Find matching bin: gene start falls within [bin_start, bin_end)
    matching_idx <- which(
      col_info$chr == gene_chr &
      col_info$start <= gene_start &
      col_info$end > gene_start
    )

    if (length(matching_idx) == 0) {
      genes_not_in_bins <- c(genes_not_in_bins, gene)
      next
    }

    # Get the column index
    col_idx <- col_info$col_idx[matching_idx[1]]

    # Add gene to this position (may have multiple genes per bin)
    if (is.null(gene_positions[[as.character(col_idx)]])) {
      gene_positions[[as.character(col_idx)]] <- gene
    } else {
      gene_positions[[as.character(col_idx)]] <- c(
        gene_positions[[as.character(col_idx)]],
        gene
      )
    }
  }

  # Warn about genes not found

if (length(genes_not_found) > 0) {
    warning("Gene(s) not found in gene_locations: ",
            paste(genes_not_found, collapse = ", "))
  }

  if (length(genes_not_in_bins) > 0) {
    warning("Gene(s) not in displayed bins (may be in filtered regions): ",
            paste(genes_not_in_bins, collapse = ", "))
  }

  # Return NULL if no genes were successfully mapped
  if (length(gene_positions) == 0) {
    return(NULL)
  }

  # Convert to at/labels format, collapsing multiple genes per bin
  at <- as.integer(names(gene_positions))
  labels <- sapply(gene_positions, function(genes) {
    paste(genes, collapse = gene_label_sep)
  })

  # Sort by position
  ord <- order(at)
  at <- at[ord]
  labels <- unname(labels[ord])

  return(list(at = at, labels = labels))
}

#' Create gene annotation for heatmap
#'
#' Creates a ComplexHeatmap annotation showing gene labels at specific genomic positions.
#'
#' @param gene_positions List with `at` (column indices) and `labels` (gene names)
#'   as returned by \code{\link{find_gene_bin_positions}}.
#' @param fontsize Font size for gene labels. Default is 10.
#' @param link_height Height of link lines in mm. Default is 5.
#'
#' @return A ComplexHeatmap HeatmapAnnotation object, or NULL if gene_positions is NULL.
#'
#' @keywords internal
make_gene_annotation <- function(gene_positions,
                                 fontsize = 10,
                                 link_height = 5) {
  if (is.null(gene_positions)) {
    return(NULL)
  }

  gene_annot <- ComplexHeatmap::HeatmapAnnotation(
    gene_labels = ComplexHeatmap::anno_mark(
      at = gene_positions$at,
      labels = gene_positions$labels,
      side = "top",
      labels_rot = 45,
      labels_gp = grid::gpar(fontsize = fontsize, fontface = "italic"),
      link_height = grid::unit(link_height, "mm"),
      padding = grid::unit(0.5, "mm"),
      extend = 0.01
    ),
    which = "column",
    show_annotation_name = FALSE
  )

  return(gene_annot)
}

make_top_annotation_gain <- function(copynumber,
                                     plotcol = "state",
                                     plotfrequency = FALSE,
                                     cutoff = NULL,
                                     maxf = NULL,
                                     frequency_height = 1.4,
                                     sv_height = 0.7,
                                     annofontsize = 10,
                                     frequency_bar_width = 0.5,
                                     SV = NULL) {
  ncells <- nrow(copynumber)

  if ((plotcol == "state" | plotcol == "copy" | plotcol == "A" | plotcol == "B") & plotfrequency == TRUE) {
    copynumbermat <- copynumber
    copynumbermat[copynumbermat == "11+"] <- "11"
    copynumbermat <- sapply(copynumbermat, as.numeric)
    f1 <- colSums(copynumbermat > cutoff, na.rm = TRUE) / ncells
    f2 <- -colSums((copynumbermat < cutoff) & (copynumbermat >= 0), na.rm = TRUE) / ncells
    if (is.null(maxf)) {
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01) {
        maxf <- 0.01
      }
    }
    ha2 <- ComplexHeatmap::columnAnnotation(
      dist2 = ComplexHeatmap::anno_barplot(
        f1,
        bar_width = frequency_bar_width,
        gp = grid::gpar(col = "#E34A33", fill = "#E34A33", lwd = 0),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c("", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2)
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        f2,
        bar_width = frequency_bar_width,
        gp = grid::gpar(col = "#3182BD", fill = "#3182BD", lwd = 0),
        axis_param = list(
          at = c(0.0, -round(maxf / 2, 2), -maxf),
          labels = c("0", "", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2)
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(frequency_height, "cm")
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
        bar_width = frequency_bar_width,
        gp = grid::gpar(
          col = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"]),
          fill = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"]),
          lwd = 0),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c("", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2)
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        matrix(data = c(f2a, f2b), ncol = 2),
        bar_width = frequency_bar_width,
        gp = grid::gpar(
          col = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"]),
          fill = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"]),
          lwd = 0),
        axis_param = list(
          at = c(0, -round(maxf / 2, 2), -maxf),
          labels = c("0", "", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2)
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(frequency_height, "cm")
    )
  }
  else if ((plotcol == "state_BAF" | plotcol == "BAF") & plotfrequency == TRUE) {
    f1 <- colSums((copynumber < 0.5) & (copynumber >= 0), na.rm = TRUE) / ncells
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
        bar_width = frequency_bar_width,
        gp = grid::gpar(col = scCNphase_colors["A-Hom"], 
                        fill = scCNphase_colors["A-Hom"],
                        lwd = 0),
        axis_param = list(
          at = c(round(maxf / 2, 2), maxf),
          labels = c("", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2, lwd = 0.3)
        ),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 = ComplexHeatmap::anno_barplot(
        f2,
        bar_width = frequency_bar_width,
        gp = grid::gpar(col = scCNphase_colors["B-Hom"], 
                        fill = scCNphase_colors["B-Hom"],
                        lwd = 0),
        axis_param = list(
          at = c(0.0, -round(maxf / 2, 2), -maxf),
          labels = c("0", "", paste0(maxf)),
          gp = grid::gpar(fontsize = annofontsize-2)
        ),
        ylim = c(-maxf, 0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(frequency_height, "cm")
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
      height = grid::unit(sv_height, "cm")
    )
  }

  return(ha2)
}

make_summary_annotations <- function(copynumber,
                                     plotcol = "state",
                                     plotmean = FALSE,
                                     plotdiversity = FALSE,
                                     mean_height = 0.7,
                                     diversity_height = 0.7,
                                     annofontsize = 10) {
  if (!plotmean && !plotdiversity) {
    return(NULL)
  }

  # Convert matrix to numeric
  mat <- copynumber
  mat[mat == "11+"] <- "11"
  mat <- suppressWarnings(apply(mat, 2, as.numeric))

  # Check if conversion produced mostly NAs (non-numeric plotcol)
  na_frac <- sum(is.na(mat)) / length(mat)
  orig_na_frac <- sum(is.na(copynumber)) / length(copynumber)
  if ((na_frac - orig_na_frac) > 0.5) {
    warning(
      "Column '", plotcol,
      "' is not numeric. Skipping mean/diversity tracks."
    )
    return(NULL)
  }

  anno_args <- list(show_annotation_name = FALSE)
  heights <- c()

  if (plotmean) {
    mean_vals <- colMeans(mat, na.rm = TRUE)
    anno_args[["mean_cn"]] <- ComplexHeatmap::anno_lines(
      mean_vals,
      gp = grid::gpar(col = "black", lwd = 1),
      add_points = FALSE,
      axis_param = list(
        side = "right",
        gp = grid::gpar(fontsize = annofontsize - 2)
      ),
      border = FALSE
    )
    heights <- c(heights, mean_height)
  }

  if (plotdiversity) {
    sd_vals <- apply(mat, 2, sd, na.rm = TRUE)
    anno_args[["diversity_cn"]] <- ComplexHeatmap::anno_lines(
      sd_vals,
      gp = grid::gpar(col = "#666666", lwd = 1),
      add_points = FALSE,
      axis_param = list(
        side = "right",
        gp = grid::gpar(fontsize = annofontsize - 2)
      ),
      border = FALSE
    )
    heights <- c(heights, diversity_height)
  }

  anno_args[["height"]] <- grid::unit(sum(heights), "cm")

  do.call(ComplexHeatmap::columnAnnotation, anno_args)
}

make_copynumber_heatmap <- function(copynumber,
                                    clones,
                                    annotations = NULL,
                                    colvals = cn_colours,
                                    legendname = "Copy Number",
                                    library_mapping = NULL,
                                    clone_pal = NULL,
                                    sample_label_idx = 1,
                                    cutoff = NULL,
                                    maxf = 1.0,
                                    plotcol = "state",
                                    plotfrequency = FALSE,
                                    frequency_height = 1.4,
                                    frequency_bar_width = 1.0,
                                    show_legend = TRUE,
                                    show_library_label = TRUE,
                                    show_clone_label = TRUE,
                                    show_clone_text = TRUE,
                                    chrlabels = TRUE,
                                    SV = NULL,
                                    labeladjust = -1,
                                    nticks = 4,
                                    Mb = TRUE,
                                    annotation_height = NULL,
                                    annofontsize = 14,
                                    na_col = "white",
                                    linkheight = 5,
                                    str_to_remove = NULL,
                                    anno_width = 0.4,
                                    rasterquality = 1,
                                    plotideogram = FALSE,
                                    genome = "hg19",
                                    ideogram_height = 0.3,
                                    legend_at = NULL,
                                    gene_annotations = NULL,
                                    gene_annotation_fontsize = NULL,
                                    gene_link_height = 5,
                                    gene_label_sep = "/",
                                    plotmean = FALSE,
                                    plotdiversity = FALSE,
                                    mean_height = 0.7,
                                    diversity_height = 0.7,
                                    ...) {

  if (class(colvals) == "function"){
    leg_params <- list(nrow = 3,
                       direction = "vertical",
                       labels_gp = grid::gpar(fontsize = annofontsize-1),
                       title_gp = grid::gpar(fontsize = annofontsize-1),
                       legend_gp = grid::gpar(fontsize = annofontsize-1))
  } else {
    # Use legend_at if provided (to exclude sentinel from legend), otherwise use all color names
    legend_values <- if (!is.null(legend_at)) legend_at else names(colvals)
    leg_params <- list(nrow = 3,
                       direction = "vertical",
                       at = legend_values,
                       labels_gp = grid::gpar(fontsize = annofontsize-1),
                       title_gp = grid::gpar(fontsize = annofontsize-1),
                       legend_gp = grid::gpar(fontsize = annofontsize-1))
  }

    # Determine which left annotation to use
  if (!is.null(annotations)) {
    left_annot <- make_left_annot_generic(
      annotations,
      show_legend = show_legend,
      annofontsize = annofontsize,
      anno_width = anno_width
    )
  } else {
    left_annot <- make_left_annot(copynumber,
      clones,
      library_mapping = library_mapping,
      clone_pal = clone_pal,
      show_clone_label = show_clone_label,
      show_clone_text = show_clone_text,
      idx = sample_label_idx,
      show_legend = show_legend,
      show_library_label = show_library_label,
      annofontsize = annofontsize,
      str_to_remove = str_to_remove,
      anno_width = anno_width
    )
  }

  # Build bottom annotation: chromosome labels + optional ideogram
  chrom_annot <- make_bottom_annot(copynumber,
    chrlabels = chrlabels,
    Mb = Mb,
    nticks = nticks,
    annotation_height = annotation_height,
    labeladjust = labeladjust,
    annofontsize = annofontsize,
    linkheight = linkheight
  )

  if (plotideogram) {
    ideogram_annot <- make_ideogram_annotation(copynumber, genome = genome,
                                               ideogram_height = ideogram_height)
    if (!is.null(ideogram_annot) && !is.null(chrom_annot)) {
      # Stack ideogram below chromosome labels
      bottom_annot <- c(ideogram_annot, chrom_annot)
    } else if (!is.null(ideogram_annot)) {
      bottom_annot <- ideogram_annot
    } else {
      bottom_annot <- chrom_annot
    }
  } else {
    bottom_annot <- chrom_annot
  }

  # Build top annotation: gene labels + frequency
  freq_annot <- make_top_annotation_gain(copynumber,
    cutoff = cutoff,
    maxf = maxf,
    plotfrequency = plotfrequency,
    plotcol = plotcol,
    SV = SV,
    frequency_height = frequency_height,
    frequency_bar_width = frequency_bar_width,
    annofontsize = annofontsize
  )

  # Create gene annotation if requested
  gene_annot <- NULL
  if (!is.null(gene_annotations)) {
    gene_fontsize <- if (!is.null(gene_annotation_fontsize)) {
      gene_annotation_fontsize
    } else {
      annofontsize
    }

    gene_positions <- find_gene_bin_positions(
      copynumber,
      gene_annotations,
      genome = genome,
      gene_label_sep = gene_label_sep
    )

    gene_annot <- make_gene_annotation(
      gene_positions,
      fontsize = gene_fontsize,
      link_height = gene_link_height
    )
  }

  # Build summary annotations (mean CN / diversity)
  summary_annot <- make_summary_annotations(
    copynumber,
    plotcol = plotcol,
    plotmean = plotmean,
    plotdiversity = plotdiversity,
    mean_height = mean_height,
    diversity_height = diversity_height,
    annofontsize = annofontsize
  )

  # Combine all top annotations (genes above frequency above summary)
  # Order: gene_annot (topmost) > freq_annot > summary_annot (closest to heatmap)
  annot_list <- Filter(Negate(is.null), list(gene_annot, freq_annot, summary_annot))
  if (length(annot_list) > 0) {
    top_annot <- Reduce(c, annot_list)
  } else {
    top_annot <- NULL
  }

  # Create the heatmap
  copynumber_hm <- ComplexHeatmap::Heatmap(
    name = legendname,
    as.matrix(copynumber),
    col = colvals,
    na_col = na_col,
    show_row_names = FALSE,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    left_annotation = left_annot,
    bottom_annotation = bottom_annot,
    heatmap_legend_param = leg_params,
    top_annotation = top_annot,
    raster_quality = rasterquality,
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
#' @param clusters data.frame assigning cells to clusters, needs the following columns `cell_id`, `clone_id` default = NULL
#' @param annotations Optional dataframe containing cell_id column and additional annotation columns
#' @param normalize_ploidy Normalize ploidy of all cells to 2
#' @param normalize_tree default = FALSE
#' @param branch_length scales branch lengths to this size, default = 2
#' @param spacer_cols number of empty columns between chromosomes, default = 20
#' @param plottree Binary value of whether to plot tree or not, default = TRUE
#' @param plotcol Which column to colour the heatmap by, should be one of "state", "state_BAF", "state_phase", "state_AS", "state_min", "copy", "BAF", "A", "B"
#' @param reorderclusters Reorder the cells according to cluster if no tree is specified
#' @param pctcells Minimum size of cluster in terms of perecentage of cells in umap clustering
#' @param library_mapping Named vector mapping library names to labels for legend
#' @param clone_pal pallette to colour clusters by
#' @param sample_label_idx default = 1
#' @param fillna Smooth over NA values, default = TRUE
#' @param frequencycutoff default = 2
#' @param frequency_bar_width Width of bars for frequency track, default = 0.5
#' @param maxf Max frequency when plotting the frequency track, default = NULL infers this from the data
#' @param plotfrequency Plot the frequency track of gains and losses across the genome
#' @param frequency_height height of the frequency track if using, default = 1.4
#' @param show_legend plot legend or not, boolean
#' @param show_library_label show library label or not, boolean
#' @param show_clone_label show clone label or not, boolean
#' @param umapmetric metric to use in umap dimensionality reduction if no clusters are specified
#' @param chrlabels include chromosome labels or not, boolean
#' @param labeladjust Adjustment for chromosome label position
#' @param SV sv data frame
#' @param seed seed for UMAP
#' @param nticks number of ticks in x-axis label when plotting a single chromosome
#' @param fillgenome Deprecated. Use plotallbins instead.
#' @param na_col colour of NA values
#' @param linkheight height of x-axis ticks
#' @param newlegendname overwrite default legend name
#' @param str_to_remove string to remove from cell_id's when plotting labels
#' @param maxCNcol max value for color scale when plotting raw data
#' @param anno_width width of left annotations
#' @param rasterquality default = 15
#' @param show_clone_text Show small inset labels next to clone/cluster annotation
#' @param widenarm Widen the copy number data table to include all bins
#' @param Mb Use Mb ticks when plotting single chromosome
#' @param annofontsize Font size to use for annotations, default = 10
#' @param annotation_height Height of the annotations
#' @param tree_width Width of phylogenetic tree, default = 4
#' @param ladderize ladderize the tree, default = TRUE, same as default in ggtree
#' @param plotallbins Include all genomic bins in the heatmap, with centromeric regions displayed as NA (light grey). Default is FALSE.
#' @param plotideogram Display chromosome ideogram (cytoband) annotation at the bottom of the heatmap. Requires plotallbins = TRUE. Default is FALSE.
#' @param genome Genome assembly to use for ideogram and centromere identification. Either "hg19" or "hg38". Default is "hg19".
#' @param centromere_col Color to use for centromeric regions when plotallbins = TRUE. Default is "#E8E8E8" (light grey).
#' @param ideogram_height Height of the ideogram annotation in cm. Default is 0.3.
#' @param gene_annotations Character vector of gene names to annotate at the top of the heatmap.
#'   Gene names must match the `ensembl_gene_symbol` column in `gene_locations` data. Default is NULL.
#' @param gene_annotation_fontsize Font size for gene annotation labels. Defaults to `annofontsize` if NULL.
#' @param gene_link_height Height of link lines connecting gene labels to positions, in mm. Default is 5.
#' @param gene_label_sep Separator used when multiple genes fall within the same genomic bin. Default is "/".
#' @param plotmean Show a mean copy number line track at the top of the heatmap. Only works with numeric plotcol values. Default is FALSE.
#' @param plotdiversity Show a copy number diversity (standard deviation) line track at the top of the heatmap. Only works with numeric plotcol values. Default is FALSE.
#' @param mean_height Height of the mean copy number track in cm. Default is 0.7.
#' @param diversity_height Height of the diversity track in cm. Default is 0.7.
#'
#' If clusters are set to NULL then the function will compute clusters using UMAP and HDBSCAN.
#' 
#' @examples
#' \dontrun{
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
                        annotations = NULL,
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
                        frequency_bar_width = 0.5,
                        maxf = NULL,
                        plotfrequency = FALSE,
                        frequency_height = 1.4,
                        show_legend = TRUE,
                        show_library_label = TRUE,
                        show_clone_label = TRUE,
                        show_clone_text = TRUE,
                        widenarm = FALSE,
                        umapmetric = "euclidean",
                        chrlabels = TRUE,
                        labeladjust = -5,
                        SV = NULL,
                        seed = NULL,
                        nticks = 4,
                        Mb = TRUE,
                        fillgenome = FALSE,
                        annotation_height = NULL,
                        annofontsize = 10,
                        na_col = "white",
                        linkheight = 2.5,
                        newlegendname = NULL,
                        str_to_remove = NULL,
                        maxCNcol = 11,
                        anno_width = 0.4,
                        rasterquality = 15,
                        tree_width = 4,
                        ladderize = TRUE,
                        plotallbins = FALSE,
                        plotideogram = FALSE,
                        genome = "hg19",
                        centromere_col = "#E8E8E8",
                        ideogram_height = 0.3,
                        gene_annotations = NULL,
                        gene_annotation_fontsize = NULL,
                        gene_link_height = 5,
                        gene_label_sep = "/",
                        plotmean = FALSE,
                        plotdiversity = FALSE,
                        mean_height = 0.7,
                        diversity_height = 0.7,
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

  if (!plotcol %in% c("state", "state_BAF", "state_phase", "state_AS", "state_min", "copy", "BAF", "B", "A", "other")) {
    stop(paste0("Column name - ", plotcol, " not available for plotting, please use one of state, copy, BAF, state_BAF, state_phase, state_AS, B or A"))
  }

  if (!plotcol %in% names(CNbins)) {
    stop(paste0("Column name - ", plotcol, " not in CNbins data frame..."))
  }

  if (plotcol == "state") {
    colvals <- cn_colours
    legendname <- "Copy Number"
  }

  if (plotcol == "B") {
    colvals <- cn_colours
    legendname <- "Copy Number\nAllele B"
  }

  if (plotcol == "A") {
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
  
  if (plotcol == "other") {
    colvals <- circlize::colorRamp2(c(0, maxCNcol / 2, maxCNcol), c(scCN_colors["CN0"], scCN_colors["CN3"], scCN_colors["CN11"]))
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

  # When plotallbins is TRUE, add centromere sentinel value to color scale
  # Centromeres use sentinel value (-999), spacers remain NA (white)
  # Store legend values before adding sentinel (to exclude sentinel from legend)
  legend_at <- NULL
  if (plotallbins) {
    sentinel_val <- CENTROMERE_SENTINEL
    if (inherits(colvals, "function")) {
      # For colorRamp2 functions, create a new function that handles the sentinel
      original_colvals <- colvals
      colvals <- function(x) {
        result <- rep(centromere_col, length(x))
        non_sentinel <- !is.na(x) & x != sentinel_val
        result[non_sentinel] <- original_colvals(x[non_sentinel])
        result[is.na(x)] <- NA
        result
      }
      # Preserve the breaks for the legend (without sentinel)
      attr(colvals, "breaks") <- attr(original_colvals, "breaks")
      class(colvals) <- class(original_colvals)
    } else {
      # Store original legend values before adding sentinel
      legend_at <- names(colvals)
      # Add sentinel value to color scale (for rendering only, not legend)
      colvals[[as.character(sentinel_val)]] <- centromere_col
    }
  }

  if (!is.null(newlegendname)){
    legendname <- newlegendname
  }

  ncells <- length(unique(CNbins$cell_id))
  
  if (!is.null(clusters) & !is.null(tree)) {
    cells_clusters <- unique(clusters$cell_id)
    cells_data <- unique(CNbins$cell_id)
    cells_tree <- unique(tree$tip.label)
    check_cells <- all(c(length(cells_tree),length(cells_clusters),length(cells_data)) == length(cells_tree))
    if (check_cells == FALSE){
      warning("Trees, clusters and copy number data have different numbers of cells, removing non-overlapping cells.")
      cells_to_keep <- intersect(intersect(cells_clusters, cells_data), cells_tree)
      CNbins <- dplyr::filter(CNbins, cell_id %in% cells_to_keep)
      clusters <- dplyr::filter(clusters, cell_id %in% cells_to_keep)
      cells_to_remove <- setdiff(cells_tree, cells_to_keep)
      tree <- ape::drop.tip(tree, cells_to_remove, collapse.singles = FALSE, trim.internal = FALSE)
      tree <- format_tree_labels(tree)
    }
  } 

  if (is.null(clusters) & !is.null(tree)) {
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
    clusters <- data.frame(cell_id = unique(CNbins$cell_id), clone_id = "0")
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
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal, ladderize = ladderize)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)

    clusters <- clustering_results$clustering %>%
      dplyr::select(cell_id, clone_id)
  }

  if (!is.null(clusters)) {
    cells_clusters <- unique(clusters$cell_id)
    cells_data <- unique(CNbins$cell_id)
    if (length(cells_data) != length(cells_clusters)){
      warning("Number of cells in clusters dataframe !=  number of cells in the bins data! Removing some cells")
      cells_to_keep <- intersect(cells_clusters, cells_data)
      CNbins <- dplyr::filter(CNbins, cell_id %in% cells_to_keep)
      clusters <- dplyr::filter(clusters, cell_id %in% cells_to_keep)
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

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal, ladderize = ladderize)
    tree_plot_dat <- tree_ggplot$data

    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
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

      tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal, ladderize = ladderize)
      tree_plot_dat <- tree_ggplot$data

      message("Creating tree...")
      tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
      ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
    }
    else if (reorderclusters == TRUE & is.null(tree)) {
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  message("Creating copy number heatmap...")

  # Handle deprecated fillgenome parameter
  if (fillgenome) {
    warning("fillgenome is deprecated. Use plotallbins = TRUE instead.")
    if (!plotallbins) {
      plotallbins <- TRUE
    }
  }

  # Validate plotideogram requires plotallbins
  if (plotideogram && !plotallbins) {
    stop("plotideogram = TRUE requires plotallbins = TRUE")
  }

  # Create copy number matrix
  # When plotallbins = TRUE, centromeres get sentinel value (colored by centromere_col in color scale)

  # Spacers remain NA (colored white by na_col)
  if (plotallbins) {
    copynumber <- createCNmatrix(CNbins, field = plotcol,
                                 plotallbins = TRUE, genome = genome)
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
  if (!is.null(annotations)) {
    if (!is.data.frame(annotations)) {
      warning("annotations is not a data.frame. Attempting conversion to data.frame")
      annotations <- as.data.frame(annotations)
    }
    annotation_idx <- match(ordered_cell_ids, annotations$cell_id)
    annotations <- annotations[annotation_idx, , drop = FALSE]
  }
  copynumber_hm <- make_copynumber_heatmap(copynumber,
    clones_formatted,
    annotations = annotations,
    colvals = colvals,
    legendname = legendname,
    library_mapping = library_mapping,
    clone_pal = clone_pal,
    sample_label_idx = sample_label_idx,
    cutoff = frequencycutoff,
    maxf = maxf,
    plotcol = plotcol,
    plotfrequency = plotfrequency,
    frequency_bar_width = frequency_bar_width,
    frequency_height = frequency_height,
    show_legend = show_legend,
    show_library_label = show_library_label,
    show_clone_label = show_clone_label,
    show_clone_text = show_clone_text,
    chrlabels = chrlabels,
    SV = SV,
    Mb = Mb,
    nticks = nticks,
    annotation_height = annotation_height,
    annofontsize = annofontsize,
    na_col = na_col,
    linkheight = linkheight,
    labeladjust = labeladjust,
    str_to_remove = str_to_remove,
    anno_width = anno_width,
    rasterquality = rasterquality,
    plotideogram = plotideogram,
    genome = genome,
    ideogram_height = ideogram_height,
    legend_at = legend_at,
    gene_annotations = gene_annotations,
    gene_annotation_fontsize = gene_annotation_fontsize,
    gene_link_height = gene_link_height,
    gene_label_sep = gene_label_sep,
    plotmean = plotmean,
    plotdiversity = plotdiversity,
    mean_height = mean_height,
    diversity_height = diversity_height,
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

  tree_ggplot <- make_tree_ggplot(tree, clusters, clone_pal = clone_pal, ladderize = ladderize)
  tree_plot_dat <- tree_ggplot$data

  message("Creating tree...")
  tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
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
    CNbins$chrarm <- paste0(CNbins$chr, coord_to_arm(CNbins$chr, CNbins$start))
    CNbins <- dplyr::left_join(CNbins, p$cl$clustering, by = "cell_id")
    CNbins <- dplyr::left_join(CNbins, p$prop)
  }

  CNbins <- CNbins %>%
    dplyr::mutate(state_BAF = ifelse(is.na(propA), "-1", state_BAF))

  plotcol <- "state_BAF"
  colvals <- cn_colours_bafstate
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
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal, ladderize = ladderize)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
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

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal, ladderize = ladderize)
    tree_plot_dat <- tree_ggplot$data

    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
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
