viterbiR <- function(emission, transition, observations) {
  emission <- t(emission)
  initial <- log(rep(1 / length(emission[, 1]), length(emission[, 1])))

  numStates <- nrow(transition)
  numObs <- length(observations)

  T1 <- matrix(data = 0, nrow = numStates, ncol = numObs)
  T2 <- matrix(data = 0, nrow = numStates, ncol = numObs)

  T1[, 1] <- initial + emission[, observations[1]]

  if (numObs > 1) {
    for (j in 2:numObs) {
      for (i in 1:numStates) {
        probs <- T1[, j - 1] + transition[, i] + emission[i, observations[j]]
        T1[i, j] <- max(probs)
        T2[i, j] <- which.max(probs)
      }
    }
  }

  # MLP = most likely path. Recover it by taking the best final state and then
  # following the stored backpointer chain (T2) from each decoded state.
  MLP <- numeric(numObs)
  MLP[numObs] <- which.max(T1[, numObs])

  if (numObs > 1) {
    for (i in numObs:2) {
      MLP[i - 1] <- T2[MLP[i], i]
    }
  }

  return(MLP - 1)
}

#' Viterbi with position-dependent transition matrices
#'
#' Reference implementation of [viterbi_pd()], kept for testing the compiled
#' version against something readable.
#'
#' @param emission Bins x states matrix of emission log likelihoods.
#' @param transition `K x K x M` array of log transition probabilities; rows are
#'   the from-state, columns the to-state.
#' @param tidx One 1-based slice index per transition, so `length(tidx)` is
#'   `length(observations) - 1`.
#' @param observations Bin indices.
#' @return Zero-based state indices, one per observation.
#' @keywords internal
viterbiR_pd <- function(emission, transition, tidx, observations) {
  emission <- t(emission)
  numStates <- dim(transition)[1]
  numObs <- length(observations)
  stopifnot(length(dim(transition)) == 3, dim(transition)[2] == numStates)
  if (numObs > 1) stopifnot(length(tidx) == numObs - 1)

  initial <- log(rep(1 / numStates, numStates))
  T1 <- matrix(0, nrow = numStates, ncol = numObs)
  T2 <- matrix(0, nrow = numStates, ncol = numObs)
  T1[, 1] <- initial + emission[, observations[1]]

  if (numObs > 1) {
    for (j in 2:numObs) {
      tr <- transition[, , tidx[j - 1]]
      for (i in 1:numStates) {
        probs <- T1[, j - 1] + tr[, i] + emission[i, observations[j]]
        T1[i, j] <- max(probs)
        T2[i, j] <- which.max(probs)
      }
    }
  }

  MLP <- numeric(numObs)
  MLP[numObs] <- which.max(T1[, numObs])
  if (numObs > 1) for (i in numObs:2) MLP[i - 1] <- T2[MLP[i], i]
  return(MLP - 1)
}
