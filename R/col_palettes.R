scCN_colors <- c(
  `CN0` = '#3182BD',
  `CN1` = '#9ECAE1',
  `CN2` = '#CCCCCC',
  `CN3` = '#FDCC8A',
  `CN4` = '#FC8D59',
  `CN5` = '#E34A33',
  `CN6` = '#B30000',
  `CN7` = '#980043',
  `CN8` = '#DD1C77',
  `CN9` = '#DF65B0',
  `CN10` = '#C994C7',
  `CN11` = '#D4B9DA'
)

scCNstate_colors <- c(
  `0` = '#3182BD',
  `1` = '#9ECAE1',
  `2` = '#CCCCCC',
  `3` = '#FDCC8A',
  `4` = '#FC8D59',
  `5` = '#E34A33',
  `6` = '#B30000',
  `7` = '#980043',
  `8` = '#DD1C77',
  `9` = '#DF65B0',
  `10` = '#C994C7',
  `11+` = '#D4B9DA'
)

scCNAS_colors <- c(
  `0|0` = '#3182BD',
  `1|0` = "#9ECAE1",
  `1|1` = "#CCCCCC",
  `2|0` = "#666666",
  `2|1` = "#FDCC8A",
  `3|0` = "#FEE2BC",
  `2|2` = "#FC8D59",
  `3|1` = "#FDC1A4",
  `4|0` = "#FB590E",
  `5` = "#E34A33",
  `6` = "#B30000",
  `7` = "#980043",
  `8` = "#DD1C77",
  `9` = "#DF65B0",
  `10` = "#C994C7",
  `11+` =  "#D4B9DA"
)

scCNminorallele_colors <- c(
  `0` = "#EAF2F8",
  `1` = "#5A9BCA",
  `2` = "#1D4E71",
  `3` = "#BD3182",
  `4` = "#5E1841",
  `5` = "#1C0713"
)

scCNphase_colors <- c(
  `A-Hom` = "#025767",
  `B-Hom` = "#A75200",
  `A-Gained` = "#53AFC0",
  `B-Gained` = "#FF9E41",
  `Balanced` = "#d5d5d4"
)

scBAFstate_colors <- c(
  `0` = "#025767",
  `0.1` = "#036F83",
  `0.2` = "#0E8BA3",
  `0.3` = "#2F99AC",
  `0.4` = "#53AFC0",
  `0.5` = "#d5d5d4",
  `0.6` = "#FFB36B",
  `0.7` = "#FF9E41",
  `0.8` = "#FF8511",
  `0.9` = "#D56800",
  `1` = "#A75200"
)

SV_colors <- c(
  `Inversion` = "darkorange3",
  `Foldback` = "#fed049",
  `Unbalanced` =  "#536162",
  `Duplication` = "#e40017",
  `Deletion` = "#78c4d4",
  `Balanced` = "#dddddd",
  `Translocation` = "#313A3A"
)

#' @export
scCN_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (scCN_colors)

  scCN_colors[cols]
}

#' @export
scCNminorallele_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (scCNminorallele_colors)

  scCNminorallele_colors[cols]
}

#' @export
scBAFstate_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (scBAFstate_colors)

  scBAFstate_colors[cols]
}

#' @export
scCNphase_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (scCNphase_colors)

  scCNphase_colors[cols]
}

#' @export
scCNstate_cols <- function(...) {
  cols <- c(...)

  if (is.null(cols))
    return (scCNstate_colors)

  scCNstate_colors[cols]
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x))), maxColorValue=255)
  col
}

scCN_palettes <- list(
  `main`= scCN_cols("CN0", "CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9", "CN10"),
  `withoutCN0` = scCN_cols("CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN7", "CN8", "CN9", "CN10"),
  `lighter1` = sapply(scCN_colors, lighten),
  `lighter2` = sapply(scCN_colors, function(x) lighten(x, factor = 1.2)),
  `darker1` = sapply(scCN_colors, function(x) lighten(x, factor = 0.9))
)

scCN_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- scCN_palettes[[palette]]

  if (reverse) pal <- rev(pal)

  colorRampPalette(pal, ...)
}

#'Colour scale constructor for scCN (single cell Copy Number) colors
#'
#' @param palette Character name of palette in scCN_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_scCN <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- scCN_pal(palette = palette, reverse = reverse)

  if (discrete) {
    ggplot2::discrete_scale("colour", paste0("scCN_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for scCN (single cell Copy Number) colors
#'
#' @param palette Character name of palette in scCN_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_scCN <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- scCN_pal(palette = palette, reverse = reverse)

  if (discrete) {
    ggplot2::discrete_scale("fill", paste0("scCN_", palette), palette = pal, ...)
  } else {
    ggplot2::scale_fill_gradientn(colours = pal(256), ...)
  }
}
