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
      dplyr::mutate(chrarm = paste0("Chr ", chrarm)) %>%
      ggplot2::ggplot(ggplot2::aes(x = BAF)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.5) +
      ggplot2::facet_wrap(~forcats::fct_reorder(chrarm, idx), scales = "free_y") +
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
      dplyr::mutate(chr = paste0("Chr ", chr)) %>%
      ggplot2::ggplot(ggplot2::aes(x = BAF)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.5) +
      ggplot2::facet_wrap(~forcats::fct_reorder(chr, idx), scales = "free_y") +
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
