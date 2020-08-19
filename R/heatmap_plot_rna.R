#' @export
plotHeatmapBAF <- function(df, removelowqcells = TRUE){
  baf <- per_arm_baf_mat(df)

  if (removelowqcells) {
    keep <- rowSums(is.na(baf$bafperchrmat)) < floor(dim(baf$bafperchrmat)[2]/2)
    baf$bafperchrmat <- baf$bafperchrmat[keep, ]
    row.names(baf$bafperchrmat) <- names(keep)[keep]
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
