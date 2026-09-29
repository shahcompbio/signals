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
