#' @export
missegregations <- function(cn, perarm = FALSE, cutoff = 0.9){

  cn <- as.data.table(cn)

  if (perarm == FALSE){
    globalmodevalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, Mode), by=list(chr, start, end),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, ploidy := mean(state)]

    percellchrvalues <- merge(cn, globalmodevalues, by = c("chr", "start", "end"), all = TRUE, suffixes = c(".cell", ".global")) %>%
      .[, ploidy.cell := mean(state.cell), by = "cell_id"] %>%
      .[, c("state_diff",
             "state_phase_diff",
             "state_min_diff") :=
          .(state.cell != state.global,
               state_phase.cell != state_phase.global,
               state_min.cell != state_min.global)] %>%
      .[, list(state_diff_count = sum(state_diff),
               state_phase_diff_count = sum(state_phase_diff),
               state_min_diff_count = sum(state_min_diff),
               state.cell = median(state.cell),
               state.global = median(state.global),
               state_phase.cell = Mode(state_phase.cell),
               state_phase.global = Mode(state_phase.global),
               state_min.cell = Mode(state_min.cell),
               state_min.global = Mode(state_min.global),
               nbins = .N), by = c("chr", "cell_id", "ploidy.cell", "ploidy")] %>%
      .[, c("state_diff_freq",
            "state_phase_diff_freq",
            "state_min_diff_freq") :=
            list(state_diff_count / nbins,
            state_phase_diff_count / nbins,
            state_min_diff_count / nbins)] %>%
      .[, c("state",
            "phase",
            "state_min") :=
          list(fifelse(state_diff_freq > cutoff, 1, 0),
               fifelse(state_phase_diff_freq > cutoff, 1, 0),
               fifelse(state_min_diff_freq > cutoff, 1, 0))]

  } else{
    cn$arm <- coord_to_arm(cn$chr, cn$start)

    globalmodevalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, Mode), by=list(chr, arm, start, end),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, ploidy := mean(state)]

    percellchrvalues <- merge(cn, globalmodevalues, by = c("chr", "arm", "start", "end"), all = TRUE, suffixes = c(".cell", ".global")) %>%
      .[, ploidy.cell := mean(state.cell), by = "cell_id"] %>%
      .[, c("state_diff",
            "state_phase_diff",
            "state_min_diff") :=
          .(state.cell != state.global,
            state_phase.cell != state_phase.global,
            state_min.cell != state_min.global)] %>%
      .[, list(state_diff_count = sum(state_diff),
               state_phase_diff_count = sum(state_phase_diff),
               state_min_diff_count = sum(state_min_diff),
               state.cell = median(state.cell),
               state.global = median(state.global),
               state_phase.cell = Mode(state_phase.cell),
               state_phase.global = Mode(state_phase.global),
               state_min.cell = Mode(state_min.cell),
               state_min.global = Mode(state_min.global),
               nbins = .N), by = c("chr","arm", "cell_id", "ploidy.cell", "ploidy")] %>%
      .[, c("state_diff_freq",
            "state_phase_diff_freq",
            "state_min_diff_freq") :=
          list(state_diff_count / nbins,
               state_phase_diff_count / nbins,
               state_min_diff_count / nbins)] %>%
      .[, c("state",
            "phase",
            "state_min") :=
          list(fifelse(state_diff_freq > cutoff, 1, 0),
               fifelse(state_phase_diff_freq > cutoff, 1, 0),
               fifelse(state_min_diff_freq > cutoff, 1, 0))] %>%
      .[, chrarm := paste0(chr, arm)]
  }

  return(as.data.frame(percellchrvalues))
}

#' @export
missegregations_vscell <- function(cn, cellid = NULL, perarm = FALSE, cutoff = 0.9){

  if (is.null(cellid)){
    cellid <- unique(cn$cell_id)[1]
  }

  cn <- as.data.table(cn)

  if (perarm == FALSE){
    globalmodevalues <- cn %>%
      .[cell_id == cellid] %>%
      .[, ploidy := mean(state, na.rm = TRUE)]
    globalmodevalues <- subset(globalmodevalues, select = c("chr", "start", "end", "state", "state_phase", "state_AS_phased", "state_min", "ploidy"))

    percellchrvalues <- merge(cn, globalmodevalues, by = c("chr", "start", "end"), all = TRUE, suffixes = c(".cell", ".global")) %>%
      .[, ploidy.cell := mean(state.cell), by = "cell_id"] %>%
      .[, c("state_diff",
            "state_phase_diff",
            "state_min_diff") :=
          .(state.cell != state.global,
            state_phase.cell != state_phase.global,
            state_min.cell != state_min.global)] %>%
      .[, ploidy := mean(ploidy, na.rm = TRUE)] %>%
      .[, list(state_diff_count = sum(state_diff, na.rm = TRUE),
               state_phase_diff_count = sum(state_phase_diff, na.rm = TRUE),
               state_min_diff_count = sum(state_min_diff, na.rm = TRUE),
               state.cell = median(state.cell, na.rm = TRUE),
               state.global = median(state.global, na.rm = TRUE),
               state_phase.cell = Mode(state_phase.cell),
               state_phase.global = Mode(state_phase.global),
               state_min.cell = Mode(state_min.cell),
               state_min.global = Mode(state_min.global),
               nbins = .N), by = c("chr", "cell_id", "ploidy.cell", "ploidy")] %>%
      .[, c("state_diff_freq",
            "state_phase_diff_freq",
            "state_min_diff_freq") :=
          list(state_diff_count / nbins,
               state_phase_diff_count / nbins,
               state_min_diff_count / nbins)] %>%
      .[, c("state",
            "phase",
            "state_min") :=
          list(fifelse(state_diff_freq > cutoff, 1, 0),
               fifelse(state_phase_diff_freq > cutoff, 1, 0),
               fifelse(state_min_diff_freq > cutoff, 1, 0))]

  } else{
    cn$arm <- coord_to_arm(cn$chr, cn$start)

    globalmodevalues <- cn %>%
      .[cell_id == cellid] %>%
      .[, ploidy := mean(state, na.rm = TRUE)]
    globalmodevalues <- subset(globalmodevalues, select = c("chr", "start", "end", "arm", "state", "state_phase", "state_AS_phased", "state_min", "ploidy"))

    percellchrvalues <- merge(cn, globalmodevalues, by = c("chr", "arm", "start", "end"), all = TRUE, suffixes = c(".cell", ".global")) %>%
      .[, ploidy.cell := mean(state.cell), by = "cell_id"] %>%
      .[, c("state_diff",
            "state_phase_diff",
            "state_min_diff") :=
          .(state.cell != state.global,
            state_phase.cell != state_phase.global,
            state_min.cell != state_min.global)] %>%
      .[, ploidy := mean(ploidy, na.rm = TRUE)] %>%
      .[, list(state_diff_count = sum(state_diff, na.rm = TRUE),
               state_phase_diff_count = sum(state_phase_diff, na.rm = TRUE),
               state_min_diff_count = sum(state_min_diff, na.rm = TRUE),
               state.cell = median(state.cell, na.rm = TRUE),
               state.global = median(state.global, na.rm = TRUE),
               state_phase.cell = Mode(state_phase.cell),
               state_phase.global = Mode(state_phase.global),
               state_min.cell = Mode(state_min.cell),
               state_min.global = Mode(state_min.global),
               nbins = .N), by = c("chr","arm", "cell_id", "ploidy.cell", "ploidy")] %>%
      .[, c("state_diff_freq",
            "state_phase_diff_freq",
            "state_min_diff_freq") :=
          list(state_diff_count / nbins,
               state_phase_diff_count / nbins,
               state_min_diff_count / nbins)] %>%
      .[, c("state",
            "phase",
            "state_min") :=
          list(fifelse(state_diff_freq > cutoff, 1, 0),
               fifelse(state_phase_diff_freq > cutoff, 1, 0),
               fifelse(state_min_diff_freq > cutoff, 1, 0))] %>%
      .[, chrarm := paste0(chr, arm)]
  }

  return(as.data.frame(percellchrvalues))
}


plotmissegs <- function(ms, arm = FALSE, yaxis = "Counts"){

  ncells <- length(unique(ms$cell_id))

  if (arm == TRUE){
    chridx <- data.frame(chrarm = paste0(rep(c(paste0(seq(1:22)), "X", "Y"), each = 2), rep(c("p", "q"), 24))) %>%
      dplyr::mutate(idx = 1:n())

    mscounts <- ms %>%
      dplyr::group_by(chrarm, state_diff) %>%
      dplyr::summarise(n = dplyr::n()) %>%
      tidyr::pivot_wider(names_from = "state_diff", values_from = n) %>%
      dplyr::select(-`NA`) %>%
      tidyr::replace_na(list(gain = 0, loss = 0)) %>%
      dplyr::mutate(loss = -loss) %>%
      tidyr::pivot_longer(cols = c(gain, loss), values_to = "Counts") %>%
      dplyr::mutate(Frequency = Counts / ncells)

    mscounts %>%
      dplyr::left_join(., chridx) %>%
      ggplot2::ggplot(ggplot2::aes(x=fct_reorder(chrarm, idx), fill = name)) +
      ggplot2::geom_bar(aes_string(y = yaxis), stat="identity", position="identity") +
      cowplot::theme_cowplot() +
      scale_fill_manual(values = c('#B30000', '#9ECAE1')) +
      xlab("Chromosome arm") +
      ylab(yaxis) +
      ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.title = element_blank())


  }

}
