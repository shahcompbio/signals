# cnbaf <- combineBAFCN(haplotypes, CNbins)
# ascn <- callAlleleSpecificCNHMM(cnbaf, ncores = 8)
#
# plotHeatmap(ascn, plotcol = "state_phase")

largestfrequencystate <- function(state, cutoff = 0.9){
  counts <- table(as.character(state))
  frequency <- counts / sum(counts)
  chr_state <- names(frequency)[frequency > cutoff]
  if (length(chr_state) == 0){
    chr_state <- NA
  }
  return(as.character(chr_state))
}

#' @export
missegregations <- function(cn, perarm = FALSE, cutoff = 0.9){

  if (perarm == FALSE){
    globalmodevalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, schnapps::Mode), by=list(chr, start, end),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, lapply(.SD, schnapps::Mode), by=list(chr),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")]

    percellchrvalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, function(x) {largestfrequencystate(x, cutoff = cutoff)}), by=list(chr, cell_id),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, state := as.numeric(state)]

    percellchrvalues <- percellchrvalues[globalmodevalues, on = "chr"]

    percellchrvalues[, state_diff := ifelse((state == i.state + 1) | (state == i.state - 1), state - i.state, NA)] %>%
      .[, state_diff := ifelse(state_diff == -1, "loss", "gain")] %>%
      .[, state_phase_diff := ifelse((state_phase != i.state_phase) & (!is.na(state_phase)), state_phase, NA) ] %>%
      .[, state_AS_diff := ifelse((state_AS != i.state_AS) & (!is.na(state_AS)), state_AS, NA) ] %>%
      .[, state_AS_phased_diff := ifelse((state_AS_phased != i.state_AS_phased)  & (!is.na(state_AS_phased)), state_AS_phased, NA) ]
  } else{
    cn$arm <- coord_to_arm(cn$chr, cn$start)

    globalmodevalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, schnapps::Mode), by=list(chr, arm, start, end),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, lapply(.SD, schnapps::Mode), by=list(chr, arm),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")]

    percellchrvalues <- cn %>%
      as.data.table() %>%
      .[, lapply(.SD, function(x) {largestfrequencystate(x, cutoff = cutoff)}), by=list(chr, arm, cell_id),
        .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
      .[, state := as.numeric(state)]

    percellchrvalues <- percellchrvalues[globalmodevalues, on = c("chr", "arm")]

    percellchrvalues[, state_diff := ifelse((state == i.state + 1) | (state == i.state - 1), state - i.state, NA)] %>%
      .[, state_diff := ifelse(state_diff == -1, "loss", "gain")] %>%
      .[, state_phase_diff := ifelse((state_phase != i.state_phase) & (!is.na(state_phase)), state_phase, NA) ] %>%
      .[, state_AS_diff := ifelse((state_AS != i.state_AS) & (!is.na(state_AS)), state_AS, NA) ] %>%
      .[, state_AS_phased_diff := ifelse((state_AS_phased != i.state_AS_phased)  & (!is.na(state_AS_phased)), state_AS_phased, NA) ] %>%
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

# countevents <- function(ms, field){
#   field <- paste0(state, "_diff")
#   as.data.table(ms)[, .(count = .N), by = list(chrarm, state_diff)] %>%
#     dcast(., chrarm  ~...)
# }
#
# ascn$arm <- coord_to_arm(ascn$chr, ascn$start)
#
# modevalues <- ascn %>%
#   as.data.table() %>%
#   .[, c("state") := list(schnapps::Mode(state)), list(chr, start, end)]
#
#
# globalmodevalues <- ascn %>%
#   as.data.table() %>%
#   .[, lapply(.SD, schnapps::Mode), by=list(chr, start, end),
#     .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
#   .[, lapply(.SD, schnapps::Mode), by=list(chr),
#     .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")]
#
# percellmodevalues <- ascn %>%
#   as.data.table() %>%
#   .[, lapply(.SD, schnapps::Mode), by=list(chr, cell_id),
#     .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")]
#
# percellchrvalues <- ascn %>%
#   as.data.table() %>%
#   .[, lapply(.SD, largestfrequencystate), by=list(chr, cell_id),
#     .SDcols=c("state", "state_phase","state_AS", "state_AS_phased", "state_min")] %>%
#   .[, state := as.numeric(state)]
#
# percellchrvalues <- percellchrvalues[globalmodevalues, on = "chr"]
#
# percellchrvalues[, state_diff := ifelse((state == i.state + 1) | (state == i.state - 1), state - i.state, NA)] %>%
#   .[, state_diff := ifelse(state_diff == -1, "gain", "loss")] %>%
#   .[, state_phase_diff := ifelse((state_phase != i.state_phase) & (!is.na(state_phase)), state_phase, NA) ] %>%
#   .[, state_AS_diff := ifelse((state_AS != i.state_AS) & (!is.na(state_AS)), state_AS, NA) ] %>%
#   .[, state_AS_phased_diff := ifelse((state_AS_phased != i.state_AS_phased)  & (!is.na(state_AS_phased)), state_AS_phased, NA) ]
#
# head(modevalues)
#
# percellchrvalues %>%
#   as.data.frame() %>%
#   #group_by(chr) %>%
#   filter(!is.na(state_diff)) %>%
#   ggplot(aes(x = state_diff, fill = chr, group = chr)) +
#   geom_bar(position = "dodge")
#
# percellchrvalues %>%
#   as.data.frame() %>%
#   #group_by(chr) %>%
#   filter(!is.na(state_phase_diff)) %>%
#   ggplot(aes(x = chr, fill = state_phase_diff, group = chr)) +
#   geom_bar(position = "dodge")
#
