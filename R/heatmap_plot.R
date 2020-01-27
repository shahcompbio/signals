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
  }
  copynumber[copynumber > 11] <- 11
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
  if (length(levels) <= 12) {
    pal <- RColorBrewer::brewer.pal(max(length(levels), 3), "Set3")
  } else if (length(levels) <= 20) {
    pal <- clone_palette_20
  } else {
    pal <- clone_palette_20
    print("WARNING: more clones than palette can accomodate!")
  }
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

make_tree_ggplot <- function(tree, clones) {
  if(!is.null(clones)) {
    clone_members <- get_clone_members(clones)
    tree <- ggtree::groupOTU(tree, clone_members)

    clone_levels <- gtools::mixedsort(unique(clones$clone_id))
    clone_pal <- make_clone_palette(clone_levels)
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
  pal <- RColorBrewer::brewer.pal(max(length(levels), 3), pal_name)
  names(pal) <- levels
  pal <- pal[levels]
  return(pal)
}

format_copynumber_values <- function(copynumber) {
  #copynumber[copynumber > 11] <- 11
  for(col in colnames(copynumber)) {
    values <- as.character(copynumber[, col])
    values[values == "11"] <- "11+"
    copynumber[, col] <- values
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

format_copynumber <- function(copynumber, ordered_cell_ids, spacer_cols=20) {
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

  copynumber <- format_copynumber_values(copynumber)
  copynumber <- space_copynumber_columns(copynumber, spacer_cols)

  return(copynumber)
}

format_clones <- function(clones, ordered_cell_ids){
  clonesdf <- dplyr::full_join(clones, data.frame(cell_id=ordered_cell_ids))
  clonesdf[is.na(clonesdf$clone_id), "clone_id"] <- "None"
  clone_counts <- clones %>% dplyr::group_by(clone_id) %>% dplyr::summarise(count=n())

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
    width=ggplot2::unit(4, "cm"),
    which="row"
  )
  tree_annot <- ComplexHeatmap::HeatmapAnnotation(
    tree=tree_annot_func, which="row", show_annotation_name=FALSE
  )

  n_cells <- sum(tree_ggplot$data$isTip)
  tree_hm <- ComplexHeatmap::Heatmap(matrix(nc=0, nr=n_cells), left_annotation=tree_annot, ...)

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

get_library_labels <- function(cell_ids) {
  labels <- sapply(strsplit(cell_ids, "-"), function(x) {
    return(x[2])
  })
  return(labels)
}

make_left_annot <- function(copynumber, clones) {
  annot_colours <- list()

  library_labels <- get_library_labels(rownames(copynumber))
  library_levels <- gtools::mixedsort(unique(library_labels))
  annot_colours$Sample <- make_discrete_palette("Set2", library_levels)
  annot_colours$Sample <- annot_colours$Sample[!is.na(annot_colours$Sample)]

  library_legend_rows <- 10

  if(!is.null(clones)) {
    clone_levels <- unique(clones$clone_label)
    clone_level_none <- clone_levels[grepl("None", clone_levels)]
    clone_levels <- gtools::mixedsort(clone_levels[!grepl("None", clone_levels)])

    clone_pal <- make_clone_palette(clone_levels)
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

    clone_legend_rows <- 10
    if(length(clone_levels) > 10) {
      clone_legend_rows <- round(sqrt(length(clone_levels) * 4))
    }

    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Clone=clones$clone_label, clone_label=clone_label_generator,
      Sample=library_labels,
      col=annot_colours, show_annotation_name=c(TRUE, FALSE, TRUE),
      which="row", annotation_width=ggplot2::unit(rep(0.4, 3), "cm"),
      annotation_legend_param=list(
        Clone=list(nrow=clone_legend_rows),
        Sample=list(nrow=library_legend_rows)
      )
    )
  } else {
    left_annot <- ComplexHeatmap::HeatmapAnnotation(
      Sample=library_labels, col=annot_colours,
      which="row", simple_anno_size=ggplot2::unit(0.4, "cm"),
      annotation_legend_param=list(
        Sample=list(nrow=library_legend_rows)
      )
    )
  }

  return(left_annot)
}

make_top_annotsnv <- function(mutgroups){

}

get_chrom_label_pos <- function(copynumber) {
  chrom_label_pos <- c()
  chroms <- sapply(strsplit(colnames(copynumber), ":"), function(x) x[[1]])
  uniq_chroms <- c(as.character(1:22), "X", "Y")
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
                     link_width = ggplot2::unit(5, "mm"), link_height = link_width,
                     link_gp = lines_gp,
                     extend = ggplot2::unit(0, "mm")) {

  which = match.arg(which)[1]

  if(!is.numeric(at)) {
    stop_wrap(paste0("`at` should be numeric ", which, " index corresponding to the matrix."))
  }

  n = length(at)
  link_gp = recycle_gp(link_gp, n)
  labels_gp = recycle_gp(labels_gp, n)
  labels2index = structure(seq_along(at), names = labels)
  at2labels = structure(labels, names = at)

  if(length(extend) == 1) extend = rep(extend, 2)
  if(length(extend) > 2) extend = extend[1:2]
  if(!inherits(extend, "unit")) extend = ggplot2::unit(extend, "npc")

  height = link_width + ComplexHeatmap::max_text_width(labels, gp = labels_gp)
  width = ggplot2::unit(1, "npc")

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
    link_height = link_height - ggplot2::unit(1, "mm")
    grid.segments(pos, ggplot2::unit(rep(1, n2), "npc"), pos, ggplot2::unit(1, "npc")-rep(link_height*(1/3), n2), default.units = "native", gp = link_gp)
    grid.segments(pos, ggplot2::unit(1, "npc")-rep(link_height*(1/3), n2), h, ggplot2::unit(1, "npc")-rep(link_height*(2/3), n2), default.units = "native", gp = link_gp)
    grid.segments(h, ggplot2::unit(1, "npc")-rep(link_height*(2/3), n2), h, ggplot2::unit(1, "npc")-rep(link_height, n2), default.units = "native", gp = link_gp)
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

make_bottom_annot <- function(copynumber) {
  chrom_label_pos <- get_chrom_label_pos(copynumber)
  bottom_annot <- ComplexHeatmap::HeatmapAnnotation(chrom_labels=anno_mark(
    at=chrom_label_pos,
    labels=names(chrom_label_pos),
    side="bottom",
    padding=0.5, extend=0.01
  ), show_annotation_name=FALSE)
  return(bottom_annot)
}

make_copynumber_heatmap <- function(copynumber,
                                    clones,
                                    colvals = cn_colours,
                                    legendname = "Copy Number", ...) {
  copynumber_hm <- ComplexHeatmap::Heatmap(
    name=legendname,
    as.matrix(copynumber),
    col=colvals,
    na_col="white",
    show_row_names=FALSE,
    cluster_rows=FALSE,
    cluster_columns=FALSE,
    show_column_names=FALSE,
    bottom_annotation=make_bottom_annot(copynumber),
    left_annotation=make_left_annot(copynumber, clones),
    heatmap_legend_param=list(nrow=4),
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
plotHeatmap <- function(CNbins,
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
                        ...){

  if (!plotcol %in% names(CNbins)){
    warning(paste0("Column name - ", plotcol, " not in CNbins data frame, using state to plot heatmap"))
    plotcol <- "state"
  }

  if (plotcol == "state"){
    colvals <- cn_colours
    legendname <- "Copy Number"
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
    clustering_results <- umap_clustering(CNbins, minPts = max(round(pctcells * ncells), 2), ...)
    tree <- clustering_results$tree
    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clustering_results$clusters))
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

    tree_ggplot <- make_tree_ggplot(tree, as.data.frame(clusters))
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

  copynumber_hm <- make_copynumber_heatmap(copynumber,
                                           format_clones(as.data.frame(clusters), ordered_cell_ids),
                                           colvals = colvals,
                                           legendname = legendname)
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
    cols <- circlize::colorRamp2(c(0, 0.0001, quantile(as.matrix(muts), 0.95)), c("#7EA5EA", "#FFFFFF", "#9A2E1C"))
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








