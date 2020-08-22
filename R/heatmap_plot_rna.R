#' @export
plotHeatmapBAF <- function(df, removelowqcells = TRUE, samplecells = NULL){
  baf <- per_arm_baf_mat(df)

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

  col_fun = circlize::colorRamp2(c(0, 0.5, 1), c(scCNphase_colors["A-LOH"], "white", scCNphase_colors["B-LOH"]))

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
