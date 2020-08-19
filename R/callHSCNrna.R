get_states_dna <- function(hscn){
  possible_states <- hscn %>%
    dplyr::filter(chr != "Y") %>%
    dplyr::mutate(arm = coord_to_arm(chr, start, mergesmallarms = TRUE)) %>%
    dplyr::mutate(chrarm = paste0(chr, arm)) %>%
    dplyr::group_by(chr, arm, chrarm, cell_id) %>%
    dplyr::summarise(state_AS_phased = Mode(state_AS_phased),
                     Maj = Mode(Maj),
                     Min = Mode(Min)) %>%
    dplyr::group_by(chr, arm, chrarm, state_AS_phased, Min, Maj) %>%
    dplyr::summarise(n = n(), f = sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(chr, arm, chrarm) %>%
    dplyr::mutate(f = n / sum(n)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(f > 0.1)

  return(possible_states)
}

#' @export
assign_states <- function(haps, hscn){

  perchrlist <- per_arm_baf_mat(haps)
  possible_states <- get_states_dna(hscn)
  perchr <- dplyr::left_join(perchrlist$bafperchr, possible_states) %>%
    dplyr::arrange(cell_id, chrarm, state_AS_phased) %>%
    as.data.table() %>%
    .[, prob := Min / (Min + Maj)] %>%
    .[, L := dbinom(size = total, x = alleleB, prob = prob)] %>%
    .[, state := Min + Maj]
  perchr <- perchr[perchr[, .I[which.max(L)], by=.(chrarm, cell_id)]$V1] %>%
    .[, state_AS_phased := paste0(Maj, "|", Min)] %>%
    .[, state_AS := paste0(pmax(state - Min, Min), "|", pmin(state - Min, Min))] %>%
    .[, state_min := pmin(Maj, Min)] %>%
    .[, state_AS := ifelse(state > 4, state, state_AS)] %>%
    .[, LOH := ifelse(state_min == 0, "LOH", "NO")] %>%
    .[, phase := c("Balanced", "A", "B")[1 +
                                           1 * ((Min < Maj)) +
                                           2 * ((Min > Maj))]] %>%
    .[, state_phase := c("Balanced", "A-Gained", "B-Gained", "A-LOH", "B-LOH")[1 +
                                                                                 1 * ((Min < Maj) & (Min != 0)) +
                                                                                 2 * ((Min > Maj) & (Maj != 0)) +
                                                                                 3 * ((Min < Maj) & (Min == 0)) +
                                                                                 4 * ((Min > Maj) & (Maj == 0))]
      ] %>%
    .[order(cell_id, chr)] %>%
    .[, state_BAF := round((Min / state)/0.1)*0.1] %>%
    .[, state_BAF := fifelse(is.nan(state_BAF), 0.5, state_BAF)] %>%
    .[, start := data.table::fifelse(arm == "p", 1, 10)] %>%
    .[, end := data.table::fifelse(arm == "p", 2, 11)] %>%
    .[, copy := state]

  perchr <- as.data.frame(perchr) %>%
    dplyr::select(-L, -n, -f, -prob) %>%
    as.data.frame()

  return(perchr)
}
