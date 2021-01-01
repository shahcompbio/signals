viterbiR <- function(emission, transition, observations) {

  emission <- t(emission)
  initial <- log(rep(1 / length(emission[, 1]), length(emission[, 1])))

  numStates <- nrow(transition)
  numObs <- length(observations)

  T1 <- matrix(data=0, nrow=numStates, ncol=numObs)
  T2 <- matrix(data=0, nrow=numStates, ncol=numObs)
  firstObs <- observations[1]

  T1[, 1] = initial + emission[, observations[1]]

  for (j in 2:length(observations)){
    for (i in 1:numStates){
      probs <- T1[ , j - 1] + transition[ ,i] + emission[i, observations[j]]
      T1[i, j] <- max(probs)
      T2[i, j] <- which.max(probs)
    }
  }

  # MLP = most likely path
  MLP <- numeric(numObs)

  MLP[numObs] <- T2[which.max(T1[,numObs]), numObs]

  # backtrace
  for (i in numObs:2) {
    zm <- which.max(T1[,i])
    MLP[i-1] <- T2[zm,i]
  }

  return(MLP - 1)
}
