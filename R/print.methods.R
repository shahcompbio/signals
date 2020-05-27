#' @export
is.hscn <- function(x) inherits(x, "hscn")

#' @export
is.ascn <- function(x) inherits(x, "ascn")

#' @export
print.hscn = function(x, ...) {
  stopifnot(is.hscn(x))
  cat("Haplotype specific copy number object \n \n")
  cat(paste0("Number of cells: ", length(unique(x$data$cell_id)), "\n"))
  cat(paste0("Bin size: ", (x$data$end[1] - x$data$start[1] + 1) / 1e6, " Mb \n"))
  cat(paste0("Inferred LOH error rate: ", round(x$loherror, 3)), "\n")
  cat(paste0("Emission model for HMM: ", x$likelihood$likelihood), "\n")
  if (x$likelihood$likelihood == "betabinomial"){
    cat(paste0("\t Inferred over dispersion: ", round(x$likelihood$rho, 4), "\n"))
    cat(paste0("\t Tarones Z score: ", round(x$likelihood$taronesZ, 3), "\n"))
  }
  cat(paste0("Average distance from median to expected BAF = ", round(x$qc_summary$summary, 4), " \n"))
}

#' @export
print.ascn = function(x, ...) {
  stopifnot(is.ascn(x))
  cat("Allele specific copy number object \n \n")
  cat(paste0("Number of cells: ", length(unique(x$data$cell_id)), "\n"))
  cat(paste0("Bin size : ", (x$data$end[1] - x$data$start[1] + 1) / 1e6, " Mb \n"))
  cat(paste0("Inferred LOH error rate : ", round(x$loherror, 3)), "\n")
  cat(paste0("Emission model for HMM: ", x$likelihood$likelihood), "\n")
  if (x$likelihood$likelihood == "betabinomial"){
    cat(paste0("\t Inferred over dispersion: ", round(x$likelihood$rho, 4), "\n"))
    cat(paste0("\t Tarones Z score: ", round(x$likelihood$taronesZ, 3), "\n"))
  }
  cat(paste0("Average distance from median to expected BAF = ", round(x$qc_summary$summary, 4), " \n"))
}
