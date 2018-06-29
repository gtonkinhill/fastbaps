#include <RcppArmadillo.h>
#include <math.h>
#include "log_add_sub.hpp"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List tree_llk(List data, arma::imat merges) {

  const arma::sp_mat& snp_matrix = data[0];

  NumericVector consensus = data[1];
  NumericMatrix prior = data[2];

  int n_isolates = snp_matrix.n_cols;
  int n_snps = snp_matrix.n_rows;
  int n_merges = merges.n_rows;

  unsigned int i,j,m,partition_length;
  int count = 0;

  arma::umat partition_counts = arma::zeros<arma::umat>(n_snps, 4);
  std::vector<arma::sp_umat> sp_partition_counts;

  // Pre-compute lgamma values up to three decimal places for each possible count lgamma(x+decimal)
  int lg_size = ceil(max(prior)*1000)+2;
  arma::dmat pre_lgamma = arma::dmat(n_isolates+2, lg_size);
  for(i=0; i<(n_isolates+2); i++){
    for(j=0; j<lg_size; j++){
      if((i==0) && (j==0)) continue;
      pre_lgamma(i,j) = lgamma(i+(j*0.001));
    }
  }
  arma::umat prior_index = arma::umat(5, n_snps);
  for(i=0; i<5; i++){
    for(j=0; j<n_snps; j++){
      prior_index(i,j) = floor(prior(i,j)/0.001);
    }
  }

  //Precompute first term in mllk
  arma::dmat term1 = arma::zeros<arma::dmat>(n_isolates+2);
  double alpha_sum = 0.0;
  for (partition_length=1; partition_length<n_isolates+2; partition_length++){
    for (j=0; j < n_snps; j++){
      alpha_sum = 0.0;
      for(i=0; i<5; i++){
        alpha_sum += prior(i,j);
      }
      term1(partition_length) -= lgamma(partition_length+alpha_sum);
    }
  }

  // Count initial alleles for each isolate into sparse matrices.
  for(i=0; i<n_isolates; i++){
    partition_counts.fill(0);
    arma::sp_mat col(snp_matrix.col(i));
    for (arma::sp_mat::iterator it = col.begin(); it != col.end(); ++it) {
      partition_counts(it.row(), *it-1) += 1;
    }
    sp_partition_counts.push_back(arma::sp_umat(partition_counts));
    count+=1;
  }

  int n_nodes = n_merges + n_isolates;
  NumericVector clust_size(n_nodes);
  std::vector<double> dk(n_nodes);
  std::vector<double> ptree(n_nodes);
  std::vector<double> mllk(n_nodes);
  std::vector<double> rk(n_nodes);
  unsigned int left, right;
  int consensus_counts[n_snps];
  double pik;

  // initilise first n_isolate elements of arrays
  for(i=0; i<n_isolates; i++){
    clust_size(i) = 1;
    dk[i] = 0.0;

    std::fill(consensus_counts, consensus_counts + n_snps, 1);
    mllk[i] = term1(1);
    for (arma::sp_umat::const_iterator it = sp_partition_counts[i].begin(); it != sp_partition_counts[i].end(); ++it) {
      consensus_counts[it.row()] -= *it;
      mllk[i]  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
    }
    for(j=0; j<n_snps; j++){
      mllk[i] += pre_lgamma(consensus_counts[j], prior_index(0, j))-pre_lgamma(0, prior_index(0, j));
    }
    ptree[i] = mllk[i];
    rk[i] = std::numeric_limits<double>::infinity();
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
    sp_partition_counts.push_back(sp_partition_counts[left] + sp_partition_counts[right]);
    count+=1;

    dk[i] = log_sum_exp(pre_lgamma(clust_size(i), 0), dk[left]+dk[right]);

    // calculate the i_th mllk
    std::fill(consensus_counts, consensus_counts + n_snps, clust_size(i));
    mllk[i] = term1(clust_size(i));
    for (arma::sp_umat::const_iterator it = sp_partition_counts[i].begin(); it != sp_partition_counts[i].end(); ++it) {
      consensus_counts[it.row()] -= *it;
      mllk[i]  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
    }
    for(j=0; j<n_snps; j++){
      mllk[i] += pre_lgamma(consensus_counts[j], prior_index(0, j))-pre_lgamma(0, prior_index(0, j));
    }

    // calculate the i_th ptree
    ptree[i] = log_sum_exp(pre_lgamma(clust_size(i), 0) + mllk[i] - dk[i],
           dk[left] + dk[right] + ptree[left] + ptree[right] - dk[i]);

    // calculate the i_th rk
    rk[i] = mllk[i] - ptree[left] - ptree[right] + pre_lgamma(clust_size(i), 0) - dk[left] - dk[right];


  }

  return(List::create(Named("dk") = wrap(dk),
                      Named("mllk") = wrap(mllk),
                      Named("rk") = wrap(rk),
                      Named("ptree") = wrap(ptree),
                      Named("n") = clust_size));

}
