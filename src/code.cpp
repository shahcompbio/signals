#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logspace_addcpp (double logx, double logy)
{
  return fmax (logx, logy) + log1p (exp (-fabs (logx - logy)));
}

// [[Rcpp::export]]
NumericVector viterbi(NumericMatrix emission, NumericMatrix transition, NumericVector observations)
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

  // Take the best final state, then follow the backpointer chain (T2) from each
  // decoded state. For numObs == 1 this returns argmax(T1(_,0)) without touching T2.
  MLP (numObs - 1) = which_max(T1 (_ , numObs - 1));

  for(int i = numObs - 1; i > 0; i--){
    MLP (i - 1) = T2 ((int) MLP (i), i);
  }

  return(MLP);
}

// Viterbi with position-dependent transition matrices.
//
// `transition` is a K x K x M array of log transition probabilities (rows =
// from-state, cols = to-state, matching viterbi() above), and `tidx` gives the
// 1-based slice to use for each of the numObs-1 transitions. Passing distinct
// matrices by index rather than one per position keeps this cheap: the allele
// aware cost depends only on the change in total copy number, which takes few
// distinct values along a chromosome.
//
// [[Rcpp::export]]
NumericVector viterbi_pd(NumericMatrix emission, NumericVector transition,
                         IntegerVector tidx, NumericVector observations)
{
  emission = transpose (emission);
  int numObs = observations.length();

  IntegerVector tdim = transition.attr("dim");
  if (tdim.length() != 3) stop("transition must be a 3-dimensional array");
  int numStates = tdim[0];
  if (tdim[1] != numStates) stop("transition slices must be square");
  if (numObs > 1 && tidx.length() != numObs - 1)
    stop("tidx must have one entry per transition (numObs - 1)");

  observations = observations - 1;

  NumericVector initial (numStates, -log(numStates));
  NumericMatrix T1 (numStates, numObs);
  NumericMatrix T2 (numStates, numObs);

  T1( _ , 0) = initial + emission( _ , observations(0));

  for(int j = 1; j < numObs; j++){
    int m = tidx(j - 1) - 1;
    if (m < 0 || m >= tdim[2]) stop("tidx out of range");
    int off = m * numStates * numStates;
    for(int i = 0; i < numStates; i++){
      NumericVector probs (numStates);
      for(int k = 0; k < numStates; k++){
        // transition[k, i, m] -> from k to i
        probs(k) = T1(k, j - 1) + transition(off + k + i * numStates)
                   + emission(i, observations(j));
      }
      T1 (i , j) = max(probs);
      T2 (i , j) = which_max(probs);
    }
  }

  NumericVector MLP (numObs);
  MLP (numObs - 1) = which_max(T1 (_ , numObs - 1));
  for(int i = numObs - 1; i > 0; i--){
    MLP (i - 1) = T2 ((int) MLP (i), i);
  }

  return(MLP);
}
