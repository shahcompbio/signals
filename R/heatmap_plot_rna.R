make_arm_matrix <- function(df){
  ord <- dplyr::distinct(df, chr, arm, chrarm) %>%
    as.data.table() %>%
    .[gtools::mixedorder(chrarm)] %>%
    .[, idx := 1:.N] %>% dplyr::as_tibble()

  baf <- dplyr::left_join(df, ord)

  baf_mat <- baf %>%
    dplyr::select(cell_id, chrarm, BAF) %>%
    tidyr::pivot_wider(names_from = "chrarm", values_from = "BAF") %>%
    as.data.frame()

  row.names(baf_mat) <- baf_mat$cell_id
  baf_mat = subset(baf_mat, select = -c(cell_id))

  idx <- dplyr::distinct(baf, chr, arm, chrarm, idx) %>% dplyr::arrange(idx)
  baf_mat <- baf_mat[, idx$chrarm]

  return(list(bafperchr = baf, bafperchrmat = baf_mat))
}

#' @export
plotHeatmapBAF <- function(df, removelowqcells = TRUE, samplecells = NULL, gen_matrix = TRUE){

  if (gen_matrix){
    baf <- per_arm_baf_mat(df)
  } else{
    baf <- make_arm_matrix(df)
  }

  if (removelowqcells) {
    keep <- rowSums(is.na(baf$bafperchrmat)) < floor(dim(baf$bafperchrmat)[2]/2)
    baf$bafperchrmat <- baf$bafperchrmat[keep, ]
    row.names(baf$bafperchrmat) <- names(keep)[keep]
  }

  if (!is.null(samplecells)) {
    cells <- baf$bafperchr %>%
      dplyr::group_by(cell_id) %>%
      dplyr::summarise(tot = sum(total)) %>%
      dplyr::arrange(desc(tot))
    samplecells <- min(samplecells, dim(cells)[1])
    cells <- cells$cell_id[1:samplecells]
    baf$bafperchrmat <- baf$bafperchrmat[cells, ]
    row.names(baf$bafperchrmat) <- cells
  }

  col_fun = circlize::colorRamp2(c(0, 0.5, 1),
                                 c(scCNphase_colors["A-LOH"], "white", scCNphase_colors["B-LOH"]))

  baf_hm <- ComplexHeatmap::Heatmap(
    name="BAF",
    as.matrix(baf$bafperchrmat),
    col=col_fun,
    na_col="grey70",
    use_raster=TRUE,
    raster_quality=5,
    show_row_names=FALSE,
    clustering_distance_rows = "pearson",
    cluster_columns=FALSE)

  return(baf_hm)
}
