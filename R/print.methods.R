#' @export
is.hscn <- function(x) inherits(x, "hscn")

#' @export
is.ascn <- function(x) inherits(x, "ascn")

#' @export
print.hscn = function(x, ...) {
  stopifnot(is.hscn(x))
  cat("Haplotype specific copy number object \n \n")
  cat(paste0("Number of cells: ", length(unique(x$data$cell_id)), "\n"))
  cat(paste0("Bin size : ", (x$data$end[1] - x$data$start[1] + 1) / 1e6, " Mb \n"))
}

#' @export
print.ascn = function(x, ...) {
  stopifnot(is.ascn(x))
  cat("Allele specific copy number object \n \n")
  cat(paste0("Number of cells: ", length(unique(x$data$cell_id)), "\n"))
  cat(paste0("Bin size : ", (x$data$end[1] - x$data$start[1] + 1) / 1e6, " Mb \n"))
}
