#' Viterbi decoding (pure R reference implementation)
#'
#' @param emission Matrix of emission log-likelihoods, observations x states.
#' @param transition Matrix of transition log-probabilities, states x states.
#' @param observations Integer vector indexing the rows of `emission` to decode.
#' @param legacy If `TRUE`, reproduce the backtrace used up to signals 0.16.0.
#'   That version seeded the final position with the *predecessor* of the best
#'   final state and then took the per-column argmax of `T1` at each step
#'   instead of following the stored backpointers, so it did not return the MAP
#'   path and largely ignored the transition penalty. Retained only to reproduce
#'   results from earlier versions.
#'
#' @return Integer vector of decoded states, 0-based.
#' @keywords internal
viterbiR <- function(emission, transition, observations, legacy = FALSE) {
  emission <- t(emission)
  initial <- log(rep(1 / length(emission[, 1]), length(emission[, 1])))

  numStates <- nrow(transition)
  numObs <- length(observations)

  T1 <- matrix(data = 0, nrow = numStates, ncol = numObs)
  T2 <- matrix(data = 0, nrow = numStates, ncol = numObs)
  firstObs <- observations[1]

  T1[, 1] <- initial + emission[, observations[1]]

  for (j in 2:length(observations)) {
    for (i in 1:numStates) {
      probs <- T1[, j - 1] + transition[, i] + emission[i, observations[j]]
      T1[i, j] <- max(probs)
      T2[i, j] <- which.max(probs)
    }
  }

  # MLP = most likely path
  MLP <- numeric(numObs)

  if (legacy) {
    MLP[numObs] <- T2[which.max(T1[, numObs]), numObs]

    for (i in numObs:2) {
      zm <- which.max(T1[, i])
      MLP[i - 1] <- T2[zm, i]
    }
  } else {
    # start at the best final state, then follow the backpointers
    s <- which.max(T1[, numObs])
    MLP[numObs] <- s

    for (i in numObs:2) {
      s <- T2[s, i]
      MLP[i - 1] <- s
    }
  }

  return(MLP - 1)
}
