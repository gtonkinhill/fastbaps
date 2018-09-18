#include <cmath>
#include <iostream>
#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

double log_sum_exp(double u, double v)
{
  if(std::isinf(u) && std::isinf(v)){
    return(-std::numeric_limits<double>::infinity());
  }
  return(std::max(u, v) + log(exp(u - std::max(u, v)) + exp(v - std::max(u, v))));
}

double log_subtract(double x, double y) {
  if(y == -std::numeric_limits<double>::infinity()){
    return x;
  }
  return x + log1p(-exp(y-x));
}

double log_sum_exp_vec(NumericVector w)
{
  double total=w[0];
  for (int i = 1; i < w.size(); ++i){
    total=log_sum_exp(total, w[i]);
  }
  return(total);
}
