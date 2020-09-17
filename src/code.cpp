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

  MLP (numObs - 1) = T2 (which_max(T1 (_ , numObs - 1)) , numObs - 1);

  for(int i = numObs - 1; i > 0; i--){
    int zm = which_max( T1 ( _, i));
    MLP (i - 1) = T2 (zm, i);
  }

  return(MLP);
}
