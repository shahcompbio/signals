#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logspace_addcpp (double logx, double logy)
{
  return fmax (logx, logy) + log1p (exp (-fabs (logx - logy)));
}

// Viterbi decoding.
//
// `legacy = true` reproduces the backtrace used up to signals 0.16.0, which was
// not a Viterbi backtrace: it seeded the final position with the *predecessor*
// of the best final state, and then took the per-column argmax of T1 at each
// step instead of following the stored backpointers. Because it ignores the
// backpointers it also largely ignores the transition penalty, so it emits
// isolated single-bin state changes that the MAP path does not contain. It is
// retained only so that results from earlier versions can be reproduced.
//
// [[Rcpp::export]]
NumericVector viterbi(NumericMatrix emission, NumericMatrix transition, NumericVector observations, bool legacy = false)
{

  emission = transpose (emission);
  int numStates = transition.nrow();
  int numObs = observations.length();

  observations = observations - 1;

  NumericVector initial  (numStates, -log(numStates));

  NumericMatrix T1 (numStates, numObs);
  NumericMatrix T2 (numStates, numObs);

  T1( _ , 0) = initial + emission( _ , observations(0));

  for(int j = 1; j < numObs; j++){
    for(int i = 0; i < numStates; i++){
      NumericVector probs = T1 ( _ , j - 1) + transition ( _ , i ) + emission (i, observations(j));
      T1 (i , j) = max(probs);
      T2 (i , j) = which_max(probs);
    }
  }

  NumericVector MLP (numObs);

  if (legacy) {
    MLP (numObs - 1) = T2 (which_max(T1 (_ , numObs - 1)) , numObs - 1);

    for(int i = numObs - 1; i > 0; i--){
      int zm = which_max( T1 ( _, i));
      MLP (i - 1) = T2 (zm, i);
    }
  } else {
    // start at the best final state, then follow the backpointers
    int s = which_max(T1 (_ , numObs - 1));
    MLP (numObs - 1) = s;

    for(int i = numObs - 1; i > 0; i--){
      s = T2 (s, i);
      MLP (i - 1) = s;
    }
  }

  return(MLP);
}
