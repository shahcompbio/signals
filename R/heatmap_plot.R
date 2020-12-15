cn_colours <- structure(
  c(
    "#3182BD", "#9ECAE1", "#CCCCCC", "#FDCC8A", "#FC8D59", "#E34A33",
    "#B30000", "#980043", "#DD1C77", "#DF65B0", "#C994C7", "#D4B9DA"
  ),
  names=c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11+")
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
  if (!is.finite(state_mode)){
    state_mode <- 2
  }
  return(state_mode)
}

normalize_cell_ploidy <- function(copynumber) {
  cell_ids <- colnames(copynumber)
  cell_ids <- cell_ids[!(cell_ids %in% c("chr", "start", "end", "width"))]

  for(cell_id in cell_ids) {
    state_mode <- calc_state_mode(copynumber[[cell_id]])
    copynumber[[cell_id]] <- as.integer(ceiling(
      copynumber[[cell_id]] / (state_mode / 2)
    ))
    copynumber[[cell_id]][copynumber[[cell_id]] > 11] <- 11
  }
  return(copynumber)
}

format_tree <- function(tree, brlen) {
  locus_tips <- grep('locus', tree$tip.label, value=TRUE)
  tree <- ape::drop.tip(tree, locus_tips)

  if (!is.null(brlen)) {
    tree <- ape::compute.brlen(tree, brlen)
  }

  tree$tip.label <- gsub('cell_', '', tree$tip.label)

  return(tree)
}

cell_order_from_tree <- function(tree, clones, cells = NULL){
  #tree <- ape::ladderize(tree, right = TRUE)
  cellorder <- tree$tip.label[tree$edge[,2]]
  clones <- as.data.frame(clones)
  row.names(clones) <- clones$cell_id
  clones <- clones[cellorder,]

  if (!is.null(cells)){
    clones <- clones[clones$cell_id %in% cells,]
  }
  clones <- clones[!is.na(clones$cell_id),]
  return(clones)
}

get_clone_members <- function(clones) {
  clone_members <- list()
  for(c in unique(clones$clone_id)) {
    if(c != "None") {
      clone_members[[c]] <- clones[clones$clone_id == c, "cell_id"]
    }
  }
  return(clone_members)
}

make_clone_palette <- function(levels) {
  if (length(levels) <= 8) {
    pal <- RColorBrewer::brewer.pal(max(length(levels), 3), "Dark2")
  } else {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), "Dark2"))(length(levels))
  }
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

make_tree_ggplot <- function(tree, clones, clone_pal = NULL) {
  if(!is.null(clones)) {
    clone_members <- get_clone_members(clones)
    tree <- ggtree::groupOTU(tree, clone_members)
    clone_levels <- gtools::mixedsort(unique(clones$clone_id))
    if (is.null(clone_pal)){
      clone_pal <- make_clone_palette(clone_levels)
    }
    clone_pal[["0"]] <- clone_none_black
    tree_aes <- ggplot2::aes(x, y, colour=group)
  } else {
    tree_aes <- ggplot2::aes(x, y)
  }

  p <- ggplot2::ggplot(tree, tree_aes) +
    ggtree::geom_tree(size=0.25) +
    ggplot2::coord_cartesian(expand=FALSE) +
    ggplot2::ylim(0.5, length(tree$tip.label) + 0.5) +
    ggplot2::theme_void()

  if(!is.null(clones)) {
    p <- p + ggplot2::scale_colour_manual(values=clone_pal)
  }

  return(p)
}

make_discrete_palette <- function(pal_name, levels) {
  if (length(levels) > 8){
    pal_name <- "Set3"
  }
  pal <- colorRampPalette(RColorBrewer::brewer.pal(max(length(levels), 3), pal_name))(length(levels))
  #pal <- RColorBrewer::brewer.pal(max(length(levels), 3), pal_name)
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

format_copynumber_values <- function(copynumber, plotcol = "state") {
  #copynumber[copynumber > 11] <- 11

  if (plotcol %in% c("BAF", "copy")){
    for(col in colnames(copynumber)) {
      values <- copynumber[, col]
      copynumber[, col] <- values
    }
  } else {
    for(col in colnames(copynumber)) {
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
    data=NA, nrow=nrow(copynumber), ncol=spacer_cols
  ))
  chrom_copynumber_dfs <- list()
  for(chrom in gtools::mixedsort(unique(chroms))) {
    chrom_copynumber <- copynumber[, chroms == chrom, drop=FALSE]
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

multi.mixedorder <- function(..., na.last = TRUE, decreasing = FALSE){
  do.call(order, c(
    lapply(list(...), function(l){
      if(is.character(l)){
        factor(l, levels=gtools::mixedsort(unique(l)))
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
                              spacer_cols=20) {
  if (!("chr" %in% colnames(copynumber))) {
    message("No chr column")
    loci <- sapply(rownames(copynumber), strsplit, "_")
    copynumber$chr <- unname(sapply(loci, '[[', 1))
    copynumber$start <- as.numeric(unname(sapply(loci, '[[', 2)))
    copynumber$end <- as.numeric(unname(sapply(loci, '[[', 3)))
    copynumber$width <- (copynumber$end - copynumber$start + 1)
  }
  copynumber$chr <- gsub("chr", "", copynumber$chr)
  copynumber <- copynumber[multi.mixedorder(copynumber$chr, copynumber$start), ]

  rownames(copynumber) <- paste0(
    copynumber$chr, ":", copynumber$start, ":", copynumber$end
  )
  copynumber <- subset(copynumber, select=-c(chr, start, end, width))
  copynumber <- as.data.frame(t(copynumber))

  copynumber <- copynumber[ordered_cell_ids, ]

  copynumber <- format_copynumber_values(copynumber, plotcol = plotcol)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)

  return(copynumber)
}

format_clones <- function(clones, ordered_cell_ids){
  clonesdf <- dplyr::full_join(clones, data.frame(cell_id=ordered_cell_ids))
  clonesdf[is.na(clonesdf$clone_id), "clone_id"] <- "None"
  clone_counts <- clones %>% dplyr::group_by(clone_id) %>% dplyr::summarise(count=dplyr::n())

  clonesdf <- dplyr::left_join(clonesdf, clone_counts) %>%
    dplyr::mutate(clone_label = paste0(clone_id, " (", count, ")")) %>%
    dplyr::select(-count) %>%
    as.data.frame()

  rownames(clonesdf) <- clonesdf$cell_id
  clonesdf <- clonesdf[ordered_cell_ids, ]
  return(clonesdf)
}

make_corrupt_tree_heatmap <- function(tree_ggplot, ...) {
  tree_annot_func = ComplexHeatmap::AnnotationFunction(
    fun=function(index) {
      pushViewport(viewport(height=1))
      grid.draw(ggplot2::ggplotGrob(tree_ggplot)$grobs[[5]])
      popViewport()
    },
    var_import=list(tree_ggplot=tree_ggplot),
    width=grid::unit(4, "cm"),
    which="row"
  )
  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree=tree_annot_func, which="row", show_annotation_name=FALSE
  )

  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- ComplexHeatmap::Heatmap(matrix(ncol=0, nrow=n_cells), left_annotation=tree_annot, ...)

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
  for(clone in unique(clones$clone_id)) {
    if(!grepl("None", clone)) {
      clone_idx <- which(clones$clone_id == clone)
      clone_idx <- find_largest_contiguous_group(clone_idx)
      clone_label_pos[[as.character(clone)]] <-
        as.integer(round(mean(clone_idx)))
    }
  }
  return(clone_label_pos)
}

get_library_labels <- function(cell_ids, idx = 1) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[idx])
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
                            show_legend = TRUE) {
  annot_colours <- list()

  if (show_clone_label == FALSE & show_library_label == FALSE){
    return(NULL)
  }

  library_labels <- get_library_labels(rownames(copynumber), idx = idx)
  if (!is.null(library_mapping)){
    library_labels <- unlist(library_mapping[library_labels])
    if (!all(library_labels %in% names(library_mapping)) == FALSE){
      warning("Not all library ids present in library to name mapping, using library IDs as annotations...")
      library_labels <- get_library_labels(rownames(copynumber))
    }
  }
  library_levels <- gtools::mixedsort(unique(library_labels))
  annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  annot_colours$Sample <- annot_colours$Sample[!is.na(annot_colours$Sample)]

  library_legend_rows <- 3

  if(!is.null(clones)) {
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])
    if (is.null(clone_pal)){
      clone_pal <- make_clone_palette(clone_levels)
    }

    if(length(clone_level_none > 0)) {
      clone_pal[[clone_level_none]] <- clone_none_black
    }
    annot_colours$Clone <- clone_pal

    clone_label_generator <- function(index) {
      clone_label_pos <- get_clone_label_pos(clones)
      y_pos <- 1 - unlist(clone_label_pos) / nrow(clones)
      grid::grid.text(
        names(clone_label_pos), 0.5, y_pos,
        just=c("centre", "centre")
      )
    }

    clone_legend_rows <- 3
    if(length(clone_levels) > 3) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 2))
    }

    if (show_library_label == TRUE & show_clone_label == TRUE){
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Clone=clones$clone_label, clone_label=clone_label_generator,
        Sample=library_labels,
        col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
        which="row", annotation_width=grid::unit(rep(0.4, 3), "cm"),
        annotation_legend_param=list(
          Clone=list(nrow=clone_legend_rows, direction = "horizontal"),
          Sample=list(nrow=library_legend_rows, direction = "horizontal")
        ),
        show_legend = show_legend
      )
    } else if (show_library_label == FALSE & show_clone_label == TRUE) {
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Clone=clones$clone_label, clone_label=clone_label_generator,
        col=annot_colours, show_annotation_name=c(TRUE, FALSE),
        which="row", annotation_width=grid::unit(rep(0.4, 2), "cm"),
        annotation_legend_param=list(
          Clone=list(nrow=clone_legend_rows)
        ),
        show_legend = show_legend
      )
    } else if (show_library_label == TRUE & show_clone_label == FALSE){
      left_annot <- ComplexHeatmap::HeatmapAnnotation(
        Sample=library_labels, col=annot_colours,
        which="row", simple_anno_size=grid::unit(0.4, "cm"),
        annotation_legend_param=list(
          Sample=list(nrow=library_legend_rows)
        ),
        show_legend = show_legend
      )
      }
  } else {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Sample=library_labels, col=annot_colours,
      which="row", simple_anno_size=grid::unit(0.4, "cm"),
      annotation_legend_param=list(
        Sample=list(nrow=library_legend_rows)
      ),
      show_legend = show_legend
    )
  }

  return(left_annot)
}

make_top_annotsnv <- function(mutgroups){

}

get_chrom_label_pos <- function(copynumber) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  chromfreq <- table(chroms)
  uniq_chroms <- names(chromfreq)[chromfreq > 1]
  uniq_chroms <- uniq_chroms[stringr::str_detect(uniq_chroms, "V", negate = TRUE)]
  for(chrom in uniq_chroms) {
    chrom_idx <- which(chroms == chrom)
    chrom_label_pos[[chrom]] <- as.integer(round(mean(chrom_idx)))
  }
  return(chrom_label_pos)
}

# From ComplexHeatmap, needed for modified anno_mark
recycle_gp = function(gp, n = 1) {
  for(i in seq_along(gp)) {
    x = gp[[i]]
    gp[[i]] = c(rep(x, floor(n/length(x))), x[seq_len(n %% length(x))])
  }
  return(gp)
}

# From ComplexHeatmap, modified
anno_mark = function(at, labels, which = c("column", "row"),
                     side = ifelse(which == "column", "top", "right"),
                     lines_gp = grid::gpar(), labels_gp = grid::gpar(), padding = 0.5,
                     link_width = grid::unit(5, "mm"), link_height = link_width,
                     link_gp = lines_gp,
                     extend = grid::unit(0, "mm")) {

  which = match.arg(which)[1]

  if(!is.numeric(at)) {
    message(paste0("`at` should be numeric ", which, " index corresponding to the matrix."))
  }

  n = length(at)
  link_gp = recycle_gp(link_gp, n)
  labels_gp = recycle_gp(labels_gp, n)
  labels2index = structure(seq_along(at), names = labels)
  at2labels = structure(labels, names = at)

  if(length(extend) == 1) extend = rep(extend, 2)
  if(length(extend) > 2) extend = extend[1:2]
  if(!inherits(extend, "unit")) extend = grid::unit(extend, "npc")

  height = link_width + ComplexHeatmap::max_text_width(labels, gp = labels_gp)
  width = grid::unit(1, "npc")

  .pos = NULL
  .scale = NULL

  column_fun = function(index) {
    n = length(index)

    # adjust at and labels
    at = intersect(index, at)
    if(length(at) == 0) {
      return(NULL)
    }
    labels = at2labels[as.character(at)]

    # labels_gp = subset_gp(labels_gp, labels2index[labels])
    # link_gp = subset_gp(link_gp, labels2index[labels])

    if(is.null(.scale)) {
      .scale = c(0.5, n+0.5)
    }
    pushViewport(viewport(yscale = c(0, 1), xscale = .scale))
    if(inherits(extend, "unit")) extend = convertWidth(extend, "native", valueOnly = TRUE)
    # text_height = convertWidth(grobHeight(textGrob(labels, gp = labels_gp))*(1+padding), "native", valueOnly = TRUE)
    text_height = convertWidth(grobWidth(textGrob(labels, gp = labels_gp))*(2+padding), "native", valueOnly = TRUE)
    if(is.null(.pos)) {
      i2 = which(index %in% at)
      pos = i2 # position of rows
    } else {
      pos = .pos[which(index %in% at)]
    }
    h1 = pos - text_height*0.5
    h2 = pos + text_height*0.5
    pos_adjusted = smartAlign(h1, h2, c(.scale[1] - extend[1], .scale[2] + extend[2]))
    h = (pos_adjusted[, 1] + pos_adjusted[, 2])/2

    n2 = length(labels)
    # grid.text(labels, h, rep(max_text_width(labels, gp = labels_gp), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
    grid.text(labels, h, rep(grobHeight(textGrob(labels, gp = labels_gp)), n2), default.units = "native", gp = labels_gp, rot = 0, just = "center")
    link_height = link_height - grid::unit(1, "mm")
    grid.segments(pos, grid::unit(rep(1, n2), "npc"), pos, grid::unit(1, "npc")-rep(link_height*(1/3), n2), default.units = "native", gp = link_gp)
    grid.segments(pos, grid::unit(1, "npc")-rep(link_height*(1/3), n2), h, grid::unit(1, "npc")-rep(link_height*(2/3), n2), default.units = "native", gp = link_gp)
    grid.segments(h, grid::unit(1, "npc")-rep(link_height*(2/3), n2), h, grid::unit(1, "npc")-rep(link_height, n2), default.units = "native", gp = link_gp)
    upViewport()
  }

  fun = column_fun

  anno = ComplexHeatmap::AnnotationFunction(
    fun = fun,
    fun_name = "anno_mark",
    which = which,
    width = width,
    height = height,
    n = -1,
    var_import = list(at, labels2index, at2labels, link_gp, labels_gp, padding, .pos, .scale,
                      side, link_width, link_height, extend),
    show_name = FALSE
  )

  # anno@subset_rule$at = subset_by_intersect

  anno@subsetable = TRUE
  return(anno)
}

make_bottom_annot <- function(copynumber, chrlabels = TRUE, filterlabels = NULL) {
  if (chrlabels[1] == FALSE){
    return(NULL)
  } else if (chrlabels[1] == TRUE){
    chrom_label_pos <- get_chrom_label_pos(copynumber)
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels=anno_mark(
      at=chrom_label_pos,
      labels=names(chrom_label_pos),
      side="bottom",
      padding=0.5, extend=0.01
    ), show_annotation_name=FALSE)
  } else {
    chrom_label_pos <- get_chrom_label_pos(copynumber)
    chrom_label_pos <- chrom_label_pos[chrlabels]
    bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels=anno_mark(
      at=chrom_label_pos,
      labels=names(chrom_label_pos),
      side="bottom",
      padding=0.5, extend=0.01
    ), show_annotation_name=FALSE)
  }
  return(bottom_annot)
}

make_top_annotation_gain <- function(copynumber,
                                     plotcol = "state",
                                     plotfrequency = FALSE,
                                     cutoff = NULL,
                                     maxf = NULL){
  ncells <- nrow(copynumber)

  if ((plotcol == "state" | plotcol == "copy") & plotfrequency == TRUE){
    f1 <- colSums(copynumber > cutoff, na.rm = TRUE) / ncells
    f2 <- -colSums(copynumber < cutoff, na.rm = TRUE) / ncells
    if (is.null(maxf)){
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01){
        maxf <- 0.01
      }
    }
    ha2 = ComplexHeatmap::columnAnnotation(
      dist2 =  ComplexHeatmap::anno_barplot(
        f1,
        bar_width = 1,
        gp =  grid::gpar(col = "#E34A33", fill = "#E34A33"),
        axis_param = list(at = c(round(maxf / 2, 2), maxf),
                          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 =  ComplexHeatmap::anno_barplot(
        f2,
        bar_width = 1,
        gp =  grid::gpar(col = "#3182BD", fill = "#3182BD"),
        axis_param = list(at = c(0.0, -round(maxf / 2, 2), -maxf),
                          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(-maxf,0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm"))
  } else if (plotcol == "state_phase" & plotfrequency == TRUE) {
    f1a <- colSums(apply(copynumber, 2, function(x) grepl("A-Gained", x))) / ncells
    f1b <- colSums(apply(copynumber, 2, function(x) grepl("A-Hom", x))) / ncells
    f2a <- -colSums(apply(copynumber, 2, function(x) grepl("B-Gained", x))) / ncells
    f2b <- -colSums(apply(copynumber, 2, function(x) grepl("B-Hom", x))) / ncells
    f1 <- f1a + f1b
    f2 <- f2a + f2b
    if (is.null(maxf)){
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01){
        maxf <- 0.01
      }
    }

    ha2 = ComplexHeatmap::columnAnnotation(
      dist2 =  ComplexHeatmap::anno_barplot(
        matrix(data = c(f1a,f1b), ncol = 2),
        bar_width = 1,
        gp =  grid::gpar(col = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"]),
                         fill = c(scCNphase_colors["A-Gained"], scCNphase_colors["A-Hom"])),
        axis_param = list(at = c(round(maxf / 2, 2), maxf),
                          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 =  ComplexHeatmap::anno_barplot(
        matrix(data = c(f2a,f2b), ncol = 2),
        bar_width = 1,
        gp =  grid::gpar(col = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"]),
                         fill = c(scCNphase_colors["B-Gained"], scCNphase_colors["B-Hom"])),
        axis_param = list(at = c(0, -round(maxf / 2, 2), -maxf),
                          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(-maxf,0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm"))
  }
  else if ((plotcol == "state_BAF" | plotcol == "BAF") & plotfrequency == TRUE){
    f1 <- colSums(copynumber < 0.5, na.rm = TRUE) / ncells
    f2 <- -colSums(copynumber > 0.5, na.rm = TRUE) / ncells
    if (is.null(maxf)){
      maxf <- ceiling(max(max(f1, max(abs(f2)))) / 0.1) * 0.1
      if (maxf < 0.01){
        maxf <- 0.01
      }
    }
    ha2 = ComplexHeatmap::columnAnnotation(
      dist2 =  ComplexHeatmap::anno_barplot(
        f1,
        bar_width = 1,
        gp =  grid::gpar(col = scCNphase_colors["A-Hom"], fill = scCNphase_colors["A-Hom"]),
        axis_param = list(at = c(round(maxf / 2, 2), maxf),
                          labels = c(paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(0, maxf),
        border = FALSE,
      ),
      dist3 =  ComplexHeatmap::anno_barplot(
        f2,
        bar_width = 1,
        gp =  grid::gpar(col = scCNphase_colors["B-Hom"], fill = scCNphase_colors["B-Hom"]),
        axis_param = list(at = c(0.0, -round(maxf / 2, 2), -maxf),
                          labels = c("0", paste0(round(maxf / 2, 2)), paste0(maxf))),
        ylim = c(-maxf,0),
        border = FALSE,
      ),
      show_annotation_name = FALSE,
      height = grid::unit(1.4, "cm"))
  }
  else {
    ha2 <- NULL
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
                                    ...) {
  copynumber_hm <- ComplexHeatmap::Heatmap(
    name=legendname,
    as.matrix(copynumber),
    col=colvals,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=make_bottom_annot(copynumber, chrlabels = chrlabels),
    left_annotation=make_left_annot(copynumber, clones,
                                    library_mapping = library_mapping, clone_pal = clone_pal, show_clone_label = show_clone_label,
                                    idx = sample_label_idx,show_legend = show_legend, show_library_label = show_library_label),
    heatmap_legend_param=list(nrow=3, direction = "vertical"),
    top_annotation = make_top_annotation_gain(copynumber, cutoff = cutoff, maxf = maxf,
                                              plotfrequency = plotfrequency, plotcol = plotcol),
    use_raster=TRUE,
    raster_quality=5,
    ...
  )
  return(copynumber_hm)
}

#' @export
createSNVmatrix <- function(SNVs, allcells = NULL, field = "mutation"){
  dfmuts <- SNVs %>%
    dplyr::mutate(mutation = ifelse(alt_counts > 0, 1, 0),
           nomutation = ifelse(alt_counts == 0, 0, 1)) %>%
    #dplyr::mutate(alt_counts = paste0(alt_counts)) %>%
    dplyr::mutate(mutid = paste(chr, as.character(start), sep = "_")) %>%
    dplyr::select_("mutid", "cell_id", field) %>%
    tidyr::spread_("mutid", field, fill = 0) %>%
    as.data.frame()

  if (!is.null(allcells)){
    message("Adding blank rows to cells that have no mutations...")
    missingcells <- data.frame(cell_id = clones$cell_id[!clones$cell_id %in% dfmuts$cell_id])
    dfmuts <- dplyr::bind_rows(dfmuts, missingcells)
  }
  rownames(dfmuts) <- dfmuts$cell_id
  dfmuts <- subset(dfmuts, select = -cell_id)

  return(dfmuts)
}

#' @export
plotHeatmap <- function(cn,
                        tree = NULL,
                        clusters = NULL,
                        normalize_ploidy = FALSE,
                        normalize_tree = FALSE,
                        branch_length = 1,
                        spacer_cols=20,
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
                        ...){

  if (is.hscn(cn) | is.ascn(cn)){
    CNbins <- cn$data
  } else{
    CNbins <- cn
  }

  if (widenarm == TRUE){
    dlpbinsarm <- dlpbins %>%
      dplyr::mutate(arm = coord_to_arm(chr, start), chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::mutate(arm = ifelse(chrarm %in% unique(CNbins$chrarm), arm, "")) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      as.data.table()

    dlpbinsarm <- data.table::rbindlist(lapply(unique(CNbins$cell_id),
                                               function(i) cbind(dlpbinsarm,
                                                                 cell_id = i))) %>%
      data.table::setkey("chr", "arm", "chrarm", "start", "end")

    CNbinst <- setkey(as.data.table(CNbins %>% dplyr::select(-start, -end)), "chr", "arm", "chrarm")
    CNbins <- dlpbinsarm[CNbinst, on = c("chr", "chrarm", "arm", "cell_id")] %>%
      .[!is.na(cell_id)] %>%
      orderdf(.)
  }

  if (!plotcol %in% c("state", "state_BAF", "state_phase", "state_AS", "state_min", "copy", "BAF")){
    stop(paste0("Column name - ", plotcol, " not available for plotting, please use one of state, copy, BAF, state_BAF, state_phase, state_AS or state_min"))
  }

  if (!plotcol %in% names(CNbins)){
    stop(paste0("Column name - ", plotcol, " not in CNbins data frame..."))
  }

  if (plotcol == "state"){
    colvals <- cn_colours
    legendname <- "Copy Number"
  }

  if (plotcol == "state_BAF"){
    colvals <- cn_colours_bafstate
    legendname <- "Allelic Imbalance"
  }

  if (plotcol == "BAF"){
    colvals = circlize::colorRamp2(c(0, 0.5, 1), c(scCNphase_colors["A-Hom"], scCNphase_colors["Balanced"], scCNphase_colors["B-Hom"]))
    legendname <- "Allelic Imbalance (Raw)"
  }

  if (plotcol == "copy"){
    colvals = circlize::colorRamp2(seq(0, 11, 1), scCN_colors)
    legendname <- "Copy"
  }

  if (plotcol == "state_AS"){
    colvals <- cn_colours_loh
    legendname <- "Allele Specific Copy Number"
  }

  if (plotcol == "state_min"){
    colvals <- cn_colours_minorallele
    legendname <- "Minor Allele Copy Number"
  }

  if (plotcol == "state_phase"){
    colvals <- cn_colours_phase
    legendname <- "Allelic Imbalance"
  }

  ncells <- length(unique(CNbins$cell_id))

  if (is.null(clusters)){
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
  } else{
    ordered_cell_ids <- paste0(clusters$cell_id)
  }

  if (is.null(tree) & is.null(clusters)){
    message("No tree or cluster information provided, clustering using HDBSCAN")
    clustering_results <- umap_clustering(CNbins,
                                          minPts = max(round(pctcells * ncells), 2),
                                          field = "copy",
                                          umapmetric = umapmetric)
    tree <- clustering_results$tree
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data
    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)

    clusters <- clustering_results$clustering %>%
      dplyr::select(cell_id, clone_id)
  }

  if (!is.null(clusters)){
    if (!"clone_id" %in% names(clusters)){
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    if (reorderclusters == TRUE){
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  if (plottree == TRUE){
    if (normalize_tree == T){
      tree <- format_tree(tree, branch_length)
    }

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
    tree_plot_dat <- tree_ggplot$data

    message("Creating tree...")
    tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
    ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
  }

  if (!is.null(clusters)){
    if (!"clone_id" %in% names(clusters)){
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    else if (reorderclusters == TRUE & !is.null(tree)){
      message("Reorder clusters dataframe according to clones using tree")
      if (normalize_tree == T){
        tree <- format_tree(tree, branch_length)
      }

      tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters), clone_pal = clone_pal)
      tree_plot_dat <- tree_ggplot$data

      message("Creating tree...")
      tree_hm <- make_corrupt_tree_heatmap(tree_ggplot)
      ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)
    }
    else if (reorderclusters == TRUE & is.null(tree)){
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  message("Creating copy number heatmap...")
  copynumber <- createCNmatrix(CNbins, field = plotcol, fillna = fillna)
  if (normalize_ploidy == T){
    message("Normalizing ploidy for each cell to 2")
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  copynumber <- format_copynumber(copynumber,
                                  ordered_cell_ids,
                                  spacer_cols = spacer_cols,
                                  plotcol = plotcol)
  clones_formatted <- format_clones(as.data.frame(clusters), ordered_cell_ids)
  if (!is.null(clone_pal)){
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
                                           ...)
  if (plottree == TRUE){
    h <- tree_hm + copynumber_hm
  } else {
    h <- copynumber_hm
  }

  return(h)
}

#' @export
plotSNVHeatmap <- function(SNVs,
                           k = 5,
                           tree = NULL,
                           clusters = NULL,
                           field = "mutation",
                           reorderclusters = FALSE){
  muts <- createSNVmatrix(SNVs, field = field)
  h <- hclust(dist(t(muts)))
  x <- cutree(h, k = k)
  mutgroups <- data.frame(MutationGroup = paste0(x), dummy = "x")
  row.names(mutgroups) <- names(x)

  muts <- muts[, h$order]
  mutgroups <- mutgroups[names(muts), ]
  mutgroups <- subset(mutgroups, select = -c(dummy))

  if (!is.null(clusters)){
    if (!"clone_id" %in% names(clusters)){
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }

    cells <- intersect(SNVs$cell_id, clusters$cell_id)
    SNVs <- dplyr::filter(SNVs, cell_id %in% cells)
    clusters <- dplyr::filter(clusters, cell_id %in% cells)

    if (reorderclusters == TRUE){
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
    }
  }

  if (is.null(clusters)){
    ordered_cell_ids <- paste0(unique(SNVs$cell_id))
  } else{
    ordered_cell_ids <- paste0(clusters$cell_id)
  }

  muts <- muts[ordered_cell_ids, ]

  colpal <- make_clone_palette(unique(mutgroups$MutationGroup))

  if (field == "mutation"){
    cols <- snv_colours
  } else {
    #cols <- circlize::colorRamp2(c(0, 0.0001, quantile(as.matrix(muts), 0.95)), c("#7EA5EA", "#FFFFFF", "#9A2E1C"))
  }

  snv_hm <- ComplexHeatmap::Heatmap(
    name="SNVs",
    as.matrix(muts),
    col=cols,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    #bottom_annotation=make_bottom_annot(copynumber),
    left_annotation=make_left_annot(muts, format_clones(clusters, ordered_cell_ids)),
    #use_raster=TRUE,
    top_annotation = HeatmapAnnotation(df = mutgroups,
                                       col = list(MutationGroup = colpal)),
    #raster_quality=1,
    heatmap_legend_param=list(nrow=4)
  )

  return(snv_hm)
}

plotHeatmapQC <- function(cn,
                        tree = NULL,
                        clusters = NULL,
                        normalize_ploidy = FALSE,
                        normalize_tree = FALSE,
                        branch_length = 1,
                        spacer_cols=20,
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
                        ...){

  CNbins <- cn$data
  p <- cn$phasing

  if (is.null(clusters)){
    clusters <- p$cl$clustering
    tree <- p$cl$tree
  }

  if ("chrarm" %in% names(p$prop)){
    arm <-TRUE
  } else {
    arm <- FALSE
  }

  if (arm == FALSE){
    CNbins <- dplyr::left_join(CNbins, p$cl$clustering)
    CNbins <- dplyr::left_join(CNbins, p$prop)
  } else{
    CNbins$chrarm <- paste0(CNbins$chr, schnapps:::coord_to_arm(CNbins$chr, CNbins$start))
    CNbins <- dplyr::left_join(CNbins, p$cl$clustering)
    CNbins <- dplyr::left_join(CNbins, p$prop)
  }

  CNbins <- CNbins %>%
   dplyr:: mutate(state_BAF = ifelse(is.na(propA), "-1", state_BAF))

  plotcol <- "state_BAF"
  colvals <- schnapps:::cn_colours_bafstate
  colvals[["-1"]] <- "gray90"
  legendname <- "Allelic Imbalance"

  ncells <- length(unique(CNbins$cell_id))

  if (is.null(clusters)){
    ordered_cell_ids <- paste0(unique(CNbins$cell_id))
  } else{
    ordered_cell_ids <- paste0(clusters$cell_id)
  }

  if (is.null(tree) & is.null(clusters)){
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

  if (!is.null(clusters)){
    if (!"clone_id" %in% names(clusters)){
      stop("No clone_id columns in clusters dataframe, you might need to rename your clusters")
    }
    if (reorderclusters == TRUE){
      message("Reorder clusters dataframe according to clones")
      clusters <- clusters[gtools::mixedorder(clusters$clone_id), ]
      ordered_cell_ids <- paste0(clusters$cell_id)
    }
  }

  if (plottree == TRUE){
    if (normalize_tree == T){
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
  if (normalize_ploidy == T){
    message("Normalizing ploidy for each cell to 2")
    copynumber <- normalize_cell_ploidy(copynumber)
  }
  copynumber <- format_copynumber(copynumber, ordered_cell_ids, spacer_cols = spacer_cols)
  clones_formatted <- format_clones(as.data.frame(clusters), ordered_cell_ids)
  if (!is.null(clone_pal)){
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
                                           sample_label_idx = sample_label_idx)
  if (plottree == TRUE){
    h <- tree_hm + copynumber_hm
  } else {
    h <- copynumber_hm
  }

  return(h)
}








