#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <omp.h>
void omp_set_num_threads(int num_threads);
int omp_get_num_threads();
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;
using namespace std;

#include "log_add_sub.hpp"
#include <cmath>

// [[Rcpp::export]]
List bhier_parallel(List data, List partitions, NumericVector d_k, int n_cores) {
  omp_set_num_threads(n_cores);

  const arma::sp_mat& snp_matrix = data[0];

  NumericVector consensus = data[1];
  NumericMatrix prior = data[2];

  const unsigned int n_partitions = partitions.size();
  const unsigned int n_isolates = snp_matrix.n_cols;
  const int n_snps = snp_matrix.n_rows;
  unsigned int i,j,p,p1,p2;
  arma::dmat rk = arma::zeros<arma::dmat>(n_partitions, n_partitions);
  arma::dmat mllk = arma::zeros<arma::dmat>(n_partitions, n_partitions);
  arma::dmat p_tree = arma::zeros<arma::dmat>(n_partitions);
  NumericVector initial_partition_llk(n_partitions);
  std::fill(rk.begin(), rk.end(), -std::numeric_limits<double>::infinity());
  std::fill(mllk.begin(), mllk.end(), -std::numeric_limits<double>::infinity());
  double d_k_t, partition_length;

  int partition_sizes[n_partitions];

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
  arma::dmat term1 = arma::zeros<arma::dmat>(n_isolates+2);
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
    partition_sizes[p] = partition_length;
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
#pragma omp parallel shared(term1, sp_partition_counts, initial_partition_llk, pre_lgamma, prior_index, p_tree, partition_sizes) default(none)
{
  unsigned int pp, k, partition_length_temp;
  int consensus_counts_temp[n_snps];
#pragma omp for
  for(pp=0; pp<n_partitions; pp++){
    partition_length_temp = partition_sizes[pp];
    std::fill(consensus_counts_temp, consensus_counts_temp + n_snps, partition_length_temp);

    p_tree(pp) = term1(partition_length_temp);

    for (arma::sp_umat::const_iterator it = sp_partition_counts[pp].begin(); it != sp_partition_counts[pp].end(); ++it) {
      consensus_counts_temp[it.row()] -= *it;
      p_tree(pp) += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
    }
    for(k=0; k<n_snps; k++){
      p_tree(pp) += pre_lgamma(consensus_counts_temp[k], prior_index(0, k))-pre_lgamma(0, prior_index(0, k));
    }
    initial_partition_llk(pp) = p_tree(pp);
  }
}

// For each pair of partitions calculate the posterior probability of combining them.
#pragma omp parallel shared(mllk, term1, sp_partition_counts, rk, pre_lgamma, prior_index, p_tree, d_k)
{
  unsigned int pp1, pp2, k, partition_length_temp;
  int consensus_counts_temp[n_snps];
#pragma omp for
  for(pp1=0; pp1<n_partitions; pp1++){
    for(pp2=(pp1+1); pp2<n_partitions; pp2++){
      arma::sp_umat res = sp_partition_counts[pp1] + sp_partition_counts[pp2];
      partition_length_temp = partition_sizes[pp1]+partition_sizes[pp2];
      std::fill(consensus_counts_temp, consensus_counts_temp + n_snps, partition_length_temp);

      mllk(pp1,pp2) = term1(partition_length_temp);

      for (arma::sp_umat::const_iterator it = res.begin(); it != res.end(); ++it) {
        consensus_counts_temp[it.row()] -= *it;
        mllk(pp1,pp2)  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
      }

      for(k=0; k<n_snps; k++){
        mllk(pp1,pp2)  += pre_lgamma(consensus_counts_temp[k], prior_index(0, k))-pre_lgamma(0, prior_index(0, k));
      }

      mllk(pp2, pp1) = mllk(pp1,pp2);

      rk(pp1,pp2) = mllk(pp1,pp2) - p_tree(pp1) - p_tree(pp2) + (pre_lgamma(partition_length_temp, 0) - d_k(pp1) - d_k(pp2));

      rk(pp2,pp1) = rk(pp1,pp2);
    }
  }
}

// Iteratively cluster the closest two clusters by mllk
IntegerVector group_mems = -seq_len(n_partitions); // Tracks group membership
NumericMatrix edges((n_partitions-1), 2);
NumericVector heights((n_partitions-1), 0.0); //hclust height output
NumericVector rk_vec((n_partitions-1), 0.0);
arma::umat used = arma::zeros<arma::umat>(n_partitions); //already clustered
NumericVector neg_inf_llk_vec(n_partitions, -std::numeric_limits<double>::infinity());
int row_index, col_index, max_i, min_i;
double temp_max;

for(j=0; j<(n_partitions-1); j++) {
  // Find max mllk and corresponding indices
  temp_max = -std::numeric_limits<double>::infinity();
  for(p1=0; p1<n_partitions; p1++){
    for(p2=(p1+1); p2<n_partitions; p2++){
      if(temp_max < rk(p1,p2)){
        temp_max = rk(p1,p2);
        row_index = p1;
        col_index = p2;
      }
    }
  }

  if(row_index > col_index){
    max_i = row_index;
    min_i = col_index;
  } else {
    max_i = col_index;
    min_i = row_index;
  }

  edges(j,0) = group_mems[min_i];
  edges(j,1) = group_mems[max_i];

  heights[j] = mllk(row_index,col_index);
  rk_vec[j] = rk(row_index,col_index);

  group_mems[min_i] = j+1;

  // Add clusters that we have merged and update partition sizes
  sp_partition_counts[min_i] = sp_partition_counts[min_i] + sp_partition_counts[max_i];
  partition_sizes[min_i] = partition_sizes[min_i] + partition_sizes[max_i];
  d_k_t = log_sum_exp(pre_lgamma(partition_sizes[min_i], 0), d_k(min_i) + d_k(max_i));
  p_tree(min_i) = log_sum_exp(pre_lgamma(partition_sizes[min_i], 0) + mllk(min_i, max_i) - d_k_t,
         d_k(min_i) + d_k(max_i) + p_tree(min_i) + p_tree(max_i) - d_k_t);
  d_k(min_i) = d_k_t;

  // Add max_i index to list of clusters we no longer need
  used[max_i] = 1;

  // update rk matrix
  partition_sizes[max_i] = -std::numeric_limits<int>::infinity();
  for(p1=0; p1<n_partitions; p1++){
    mllk(p1, max_i) = -std::numeric_limits<double>::infinity();
    mllk(max_i, p1) = -std::numeric_limits<double>::infinity();
    rk(p1, max_i) = -std::numeric_limits<double>::infinity();
    rk(max_i, p1) = -std::numeric_limits<double>::infinity();
  }

#pragma omp parallel shared(used, term1, mllk, rk, sp_partition_counts, initial_partition_llk, pre_lgamma, prior_index,d_k, p_tree) firstprivate(min_i)
{
  unsigned int pp, k, partition_length_temp;
  int consensus_counts_temp[n_snps];
#pragma omp for
  for(pp=0; pp<n_partitions; pp++){
    if((used[pp]!=1) && (pp!=min_i)){

      arma::sp_umat res = sp_partition_counts[pp] + sp_partition_counts[min_i];
      partition_length_temp = partition_sizes[pp]+partition_sizes[min_i];
      std::fill(consensus_counts_temp, consensus_counts_temp + n_snps, partition_length_temp);

      mllk(min_i,pp) = term1(partition_length_temp);

      for (arma::sp_umat::const_iterator it = res.begin(); it != res.end(); ++it) {
        consensus_counts_temp[it.row()] -= *it;
        mllk(min_i,pp)  += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
      }

      for(k=0; k<n_snps; k++){
        mllk(min_i,pp)  += pre_lgamma(consensus_counts_temp[k], prior_index(0, k))-pre_lgamma(0, prior_index(0, k));
      }

      mllk(pp,min_i) = mllk(min_i,pp);

      rk(pp,min_i) = mllk(min_i,pp) - p_tree(min_i) - p_tree(pp) + pre_lgamma(partition_length_temp, 0) - d_k(min_i) - d_k(pp);
      rk(min_i,pp) = rk(pp,min_i);
    }

  }
}

}

return(List::create(Named("initial_partition_llk") = initial_partition_llk,
                    Named("edges") = edges,
                    Named("heights")= heights,
                    Named("rk")=rk_vec));

}
