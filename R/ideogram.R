#' @export
plotideogram <- function(cn, cellid = NULL, gene.symbols = NULL, chr = NULL, maxCN = 10, version = 91, grch = "37"){

  if (!requireNamespace("GenomeInfoDb", quietly = TRUE)) {
    stop("Package \"GenomeInfoDb\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("regioneR", quietly = TRUE)) {
    stop("Package \"regioneR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("karyoploteR", quietly = TRUE)) {
    stop("Package \"karyoploteR\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package \"biomaRt\" needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (is.hscn(cn) | is.ascn(cn)){
    data <- cn$data
  } else{
    data <- cn
  }

  if (is.null(cellid)){
    cellid = data$cell_id[1]
  }

  data_cell <- data %>% dplyr::filter(cell_id == cellid)
  data_cell$col_state <- scCN_cols(paste0("CN", data_cell$state))
  data_cell$col_ASstate <- scCNphase_cols(paste0(data_cell$state_phase))

  data <- regioneR::toGRanges(data_cell[, c("chr", "start", "end", "state", "copy",
                                            "state_AS_phased","state_phase", "BAF",
                                            "col_state", "col_ASstate")])
  GenomeInfoDb::seqlevelsStyle(data) <- "UCSC"

  if (is.null(chr)){
    chr <- "auto"
  } else{
    chr <- paste0("chr", chr)
  }

  if(is.null(gene.symbols)){
    pp <- karyoploteR::getDefaultPlotParams(plot.type = 4)
    pp$data1inmargin <- 2
    pl <- expression(kp <- karyoploteR::plotKaryotype(plot.type = 4, plot.params = pp, chromosomes = chr),
      karyoploteR::kpAxis(kp, r0=0.52, r1=1.0),
      karyoploteR::kpPoints(kp, data=data, y=data$BAF, r0=0.52, r1=1.0, col = data$col_ASstate),
      karyoploteR::kpAxis(kp, tick.pos = seq(0, maxCN, 2), r0=0, r1=0.48, ymax=maxCN, ymin=0),
      karyoploteR::kpPoints(kp, data=data, y=data$copy, r0=0.0, r1=0.48, ymin = 0, ymax = maxCN, col = data$col_state))
  } else{
    ensembl <- biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=version, GRCh = grch)
    genes <- regioneR::toGRanges(biomaRt::getBM(attributes=c('chromosome_name', 'start_position', 'end_position', 'hgnc_symbol'),
                                                filters = 'hgnc_symbol', values =gene.symbols, mart = ensembl))
    GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"

    pp <- karyoploteR::getDefaultPlotParams(plot.type = 4)
    pp$data1inmargin <- 2
    pl <- expression(kp <- karyoploteR::plotKaryotype(plot.type = 4, plot.params = pp, chromosomes = chr),
      #karyoploteR::kpAddCytobandsAsLine(kp)
      karyoploteR::kpAxis(kp, r0=0.4, r1=0.75),
      karyoploteR::kpPoints(kp, data=data, y=data$BAF, r0=0.4, r1=0.75, col = data$col_ASstate),
      karyoploteR::kpAxis(kp, tick.pos = seq(0, maxCN, 2), r0=0, r1=0.35, ymax=maxCN, ymin=0),
      karyoploteR::kpPoints(kp, data=data, y=data$copy, r0=0.0, r1=0.35, ymin = 0, ymax = maxCN, col = data$col_state),
      karyoploteR::kpPlotMarkers(kp, data=genes, labels=genes$hgnc_symbol, line.color = "#555555", marker.parts = c(0.95,0.025,0.025),  r1=1.05))

  }
  eval(pl)
  return(pl)
}

