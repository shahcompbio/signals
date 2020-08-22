#' @export
per_chr_baf <- function(haps, filtern = 1, perarm = FALSE){

  if (perarm){

    chridx <- data.frame(chrarm = paste0(rep(c(paste0(seq(1:22)), "X", "Y"), each = 2), rep(c("p", "q"), 24))) %>%
      dplyr::mutate(idx = 1:n())

    hst <- haps %>%
      dplyr::mutate(arm = coord_to_arm(chr, start)) %>%
      dplyr::mutate(chrarm = paste0(chr, arm)) %>%
      dplyr::group_by(chr, arm, chrarm, cell_id) %>%
      dplyr::summarise(alleleA = sum(alleleA, na.rm = TRUE),
                       alleleB = sum(alleleB, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(BAF = alleleB / (alleleA + alleleB),
                    total = alleleB + alleleA) %>%
      dplyr::filter(total > filtern) %>%
      dplyr::left_join(chridx, by = c("chrarm")) %>%
      dplyr::mutate(chrarmf = forcats::fct_reorder(as.factor(paste0("Chr ", chrarm)), idx)) %>%
      ggplot2::ggplot(ggplot2::aes(x = BAF)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.5) +
      ggplot2::facet_wrap(~chrarmf, scales = "free_y")+
      cowplot::theme_cowplot() +
      ggplot2::scale_x_continuous(breaks = c(0.0, 0.5, 1.0)) +
      cowplot::panel_border() +
      ggplot2::geom_vline(xintercept = 0.5, lty = 2, size = 0.5, col = "firebrick4") +
      ggplot2::theme(axis.title.y=ggplot2::element_blank(),
                     axis.text.y=ggplot2::element_blank(),
                     axis.ticks.y=ggplot2::element_blank())

  } else {
    chridx <- data.frame(chr = c(paste0(1:22), "X", "Y"), idx = seq(1:24))

    hst <- haps %>%
      dplyr::group_by(chr, cell_id) %>%
      dplyr::summarise(alleleA = sum(alleleA, na.rm = TRUE),
                       alleleB = sum(alleleB, na.rm = TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(BAF = alleleB / (alleleA + alleleB),
                    total = alleleB + alleleA) %>%
      dplyr::filter(total > filtern) %>%
      dplyr::left_join(chridx, by = "chr") %>%
      dplyr::mutate(chr = forcats::fct_reorder(as.factor(paste0("Chr ", chr)), idx)) %>%
      ggplot2::ggplot(ggplot2::aes(x = BAF)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.5) +
      ggplot2::facet_wrap(~chr, scales = "free_y") +
      cowplot::theme_cowplot() +
      ggplot2::scale_x_continuous(breaks = c(0.0, 0.5, 1.0)) +
      cowplot::panel_border() +
      ggplot2::geom_vline(xintercept = 0.5, lty = 2, size = 0.5, col = "firebrick4") +
      ggplot2::theme(axis.title.y=ggplot2::element_blank(),
            axis.text.y=ggplot2::element_blank(),
            axis.ticks.y=ggplot2::element_blank())
  }

  return(hst)
}

#' @export
plot_proportions <- function(hscn_dna_arm, hscn_rna_arm, perarm = FALSE){

  prop_rna <- hscn_rna_arm %>%
    dplyr::mutate(state_AS_phased = ifelse(state > 5, "Other", state_AS_phased)) %>%
    dplyr::group_by(chr, chrarm, state_AS_phased) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(f = n / sum(n)) %>%
    dplyr::mutate(assay = "RNA")

  prop_dna <- hscn_dna_arm %>%
    dplyr::mutate(state_AS_phased = ifelse(state > 5, "Other", state_AS_phased)) %>%
    dplyr::group_by(chr, chrarm, state_AS_phased) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::mutate(f = n / sum(n)) %>%
    dplyr::mutate(assay = "DNA")

  prop <- dplyr::bind_rows(prop_rna, prop_dna) %>%
    dplyr::mutate(ord = ifelse(assay == "DNA", 1, 0)) %>%
    dplyr::mutate(chrord = ifelse(chr == "X", 24, as.numeric(as.character(chr))))

  barplot <- prop %>%
    ggplot2::ggplot(ggplot2::aes(x = assay,
               y = f,
               fill = state_AS_phased,
               alpha = forcats::fct_reorder(assay, ord))) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::facet_wrap(~forcats::fct_reorder(chrarm, chrord), nrow = 1) +
    cowplot::theme_cowplot() +
    ggplot2::theme(panel.spacing = grid::unit(0.1, "lines"),
          legend.position = "bottom",
          legend.title = ggplot2::element_blank()) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 3),
           alpha = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::xlab("") +
    ggplot2::ylab("Proportion") +
    ggplot2::scale_alpha_discrete(range = c(0.65, 1.0)) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank(),
            axis.line.x=ggplot2::element_blank()) +
    ggplot2::scale_y_continuous(expand = c(0, 0))

  return(barplot)
}
