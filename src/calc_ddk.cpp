#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <math.h>
#include "log_add_sub.hpp"


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List calc_ddk(List data, arma::imat merges) {

  const arma::sp_mat& snp_matrix = data[0];

  NumericVector consensus = data[1];
  NumericMatrix prior = data[2];

  int n_isolates = snp_matrix.n_cols;
  int n_snps = snp_matrix.n_rows;
  int n_merges = merges.n_rows;

  unsigned int i,j,m,partition_length,left,right;
  int count = 0;

  int n_nodes = n_merges + n_isolates;
  NumericVector clust_size(n_nodes);
  std::vector<double> dk(n_nodes);

  // initilise first n_isolate elements of arrays
  for(i=0; i<n_isolates; i++){
    clust_size(i) = 1;
    dk[i] = log(1);
  }

  for(m=0; m<n_merges; m++){
    // calculate the m_th count matrix
    i = m + n_isolates;
    if(merges(m,0)<0){
      left = std::abs(merges(m,0)) - 1;
    } else {
      left = n_isolates + merges(m,0) - 1;
    }
    if(merges(m,1)<0){
      right = std::abs(merges(m,1)) - 1;
    } else {
      right = n_isolates + merges(m,1) - 1;
    }

    clust_size(i) = clust_size(left) + clust_size(right);

    dk[i] = log_sum_exp(lgamma(clust_size(i)), dk[left]+dk[right]);

  }

  return(List::create(Named("dk") = wrap(dk),
                      Named("n") = clust_size));

}
