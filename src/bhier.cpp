#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "log_add_sub.hpp"
#include <math.h>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List bhier(List data, List partitions, NumericVector d_k) {
  const arma::sp_mat& snp_matrix = data[0];

  NumericVector consensus = data[1];
  NumericMatrix prior = data[2];

  int n_partitions = partitions.size();
  int n_isolates = snp_matrix.n_cols;
  int n_snps = snp_matrix.n_rows;
  size_t i,j,p,p1,p2;
  NumericMatrix rk(n_partitions, n_partitions);
  NumericMatrix mllk(n_partitions, n_partitions);
  NumericVector p_tree(n_partitions);
  NumericVector initial_partition_llk(n_partitions);
  std::fill(rk.begin(), rk.end(), -std::numeric_limits<double>::infinity());
  std::fill(mllk.begin(), mllk.end(), -std::numeric_limits<double>::infinity());
  double d_k_t, partition_length;

  NumericVector partition_sizes(n_partitions);

  arma::umat consensus_counts = arma::zeros<arma::umat>(n_snps);
  arma::umat partition_counts = arma::zeros<arma::umat>(n_snps, 4);
  std::vector<arma::sp_umat> sp_partition_counts;

  // Pre-compute lgamma values up to three decimal places for each possible count lgamma(x+decimal)
  int lg_size = ceil(max(prior)*1000)+3;
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
      prior_index(i,j) = ceil(prior(i,j)/0.001);
    }
  }

  //Precompute first term in mllk
  NumericVector term1(n_isolates+2, 0.0);
  double alpha_sum = 0.0;
  for (partition_length=1; partition_length<n_isolates+2; partition_length++){
    for (j=0; j < n_snps; j++){
      alpha_sum = 0.0;
      for(i=0; i<5; i++){
        alpha_sum += prior(i,j);
      }
      term1(partition_length) += lgamma(alpha_sum) - lgamma(partition_length+alpha_sum);
    }
  }


  // Count initial alleles for each partition into sparse matrices.
  for(p=0; p<n_partitions; p++){
    IntegerVector partition = partitions(p);
    partition_length = partition.size();
    partition_sizes(p) = partition_length;
    partition_counts.fill(0);

    for (i=0; i<partition_length; i++){
      arma::sp_mat col(snp_matrix.col(partition(i)-1));
      for (arma::sp_mat::iterator it = col.begin(); it != col.end(); ++it) {
        partition_counts(it.row(), *it-1) += 1;
      }
    }

    sp_partition_counts.push_back(arma::sp_umat(partition_counts));
  }


  // Calculate the llk of each partition
  for(p=0; p<n_partitions; p++){
    partition_length = partition_sizes(p);
    consensus_counts.fill(partition_length);
    p_tree(p) = term1(partition_length); //-n_snps*pre_lgamma(partition_length+1, 0);

    for (arma::sp_umat::const_iterator it = sp_partition_counts[p].begin(); it != sp_partition_counts[p].end(); ++it) {
      consensus_counts[it.row()] -= *it;
      p_tree(p) += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
    }
    for(i=0; i<n_snps; i++){
      p_tree(p) += pre_lgamma(consensus_counts[i], prior_index(0, i))-pre_lgamma(0, prior_index(0, i));
    }
    initial_partition_llk(p) = p_tree(p);
  }

  // For each pair of partitions calculate the posteriror probability of combining them.
  for(p1=0; p1<n_partitions; p1++){
    for(p2=(p1+1); p2<n_partitions; p2++){
      arma::sp_umat res = sp_partition_counts[p1] + sp_partition_counts[p2];
      partition_length = partition_sizes(p1)+partition_sizes(p2);
      consensus_counts.fill(partition_length);

      mllk(p1,p2) = term1(partition_length); //-n_snps*pre_lgamma(partition_length+1, 0);

      for (arma::sp_umat::const_iterator it = res.begin(); it != res.end(); ++it) {
        consensus_counts[it.row()] -= *it;
        mllk(p1,p2)  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
      }

      for(i=0; i<n_snps; i++){
        mllk(p1,p2)  += pre_lgamma(consensus_counts[i], prior_index(0, i))-pre_lgamma(0, prior_index(0, i));
      }

      mllk(p2, p1) = mllk(p1,p2);

      rk(p1,p2) = mllk(p1,p2) - p_tree(p1) - p_tree(p2) + (pre_lgamma(partition_length, 0) - d_k(p1) - d_k(p2));

      rk(p2,p1) = rk(p1,p2);

    }
  }

  // Iteratively cluster the closest two clusters by mllk
  IntegerVector group_mems = -seq_len(n_partitions); // Tracks group membership
  NumericMatrix edges((n_partitions-1), 2);
  NumericVector heights((n_partitions-1), 0.0); //hclust height output
  NumericVector rk_vec((n_partitions-1), 0.0);
  arma::umat used = arma::zeros<arma::umat>(n_partitions); //already clustered
  NumericVector neg_inf_llk_vec(n_partitions, -std::numeric_limits<double>::infinity());
  int max_index, row_index, col_index, max_i, min_i;

  for(j=0; j<(n_partitions-1); j++) {
    // Find max mllk and corresponding indices
    max_index = which_max(rk);

    row_index = max_index % n_partitions;
    col_index = max_index / n_partitions;

    if(row_index > col_index){
      max_i = row_index;
      min_i = col_index;
    } else {
      max_i = col_index;
      min_i = row_index;
    }

    edges(j,0) = group_mems[min_i];
    edges(j,1) = group_mems[max_i];

    heights[j] = mllk[max_index];
    rk_vec[j] = rk[max_index];

    group_mems[min_i] = j+1;

    // Rcout << "min_i: " << min_i << std::endl;
    // Rcout << "max_i: " << max_i << std::endl;

    // Add clusters that we have merged and update partition sizes
    sp_partition_counts[min_i] = sp_partition_counts[min_i] + sp_partition_counts[max_i];
    partition_sizes(min_i) = partition_sizes(min_i) + partition_sizes(max_i);
    d_k_t = log_sum_exp(pre_lgamma(partition_sizes(min_i), 0), d_k(min_i) + d_k(max_i));
    p_tree(min_i) = log_sum_exp(pre_lgamma(partition_sizes(min_i), 0) + mllk(min_i, max_i) - d_k_t,
           d_k(min_i) + d_k(max_i) + p_tree(min_i) + p_tree(max_i) - d_k_t);
    d_k(min_i) = d_k_t;


    // Add max_i index to list of clusters we no longer need
    used[max_i] = 1;

    // update rk matrix
    partition_sizes(max_i) = -std::numeric_limits<double>::infinity();
    mllk(_, max_i) = neg_inf_llk_vec;
    mllk(max_i, _) = neg_inf_llk_vec;
    rk(_, max_i) = neg_inf_llk_vec;
    rk(max_i, _) = neg_inf_llk_vec;

    for(p=0; p<n_partitions; p++){
      if(used[p]==1)  continue;
      if(p==min_i) continue;

      arma::sp_umat res = sp_partition_counts[p] + sp_partition_counts[min_i];
      partition_length = partition_sizes(p)+partition_sizes(min_i);
      consensus_counts.fill(partition_length);

      mllk(min_i,p) = term1(partition_length); //-n_snps*pre_lgamma(partition_length+1, 0);;

      for (arma::sp_umat::const_iterator it = res.begin(); it != res.end(); ++it) {
        consensus_counts[it.row()] -= *it;
        mllk(min_i,p)  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
      }

      for(i=0; i<n_snps; i++){
        mllk(min_i,p)  += pre_lgamma(consensus_counts[i], prior_index(0, i))-pre_lgamma(0, prior_index(0, i));
      }

      mllk(p,min_i) = mllk(min_i,p);

      rk(p,min_i) = mllk(min_i,p) - p_tree(min_i) - p_tree(p) + pre_lgamma(partition_length, 0) - d_k(min_i) - d_k(p);
      rk(min_i,p) = rk(p,min_i);

    }
  }

  return(List::create(Named("initial_partition_llk") = initial_partition_llk,
                      Named("edges") = edges,
                      Named("heights")= heights,
                      Named("rk")=rk_vec));

}
