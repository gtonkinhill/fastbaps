#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]


#include <math.h>


using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List part_llks(List data, List partitions) {

  const arma::sp_mat& snp_matrix = data[0];

  NumericVector consensus = data[1];
  NumericMatrix prior = data[2];

  int n_partitions = partitions.size();
  int n_isolates = snp_matrix.n_cols;
  int n_snps = snp_matrix.n_rows;
  int i,j,p,a,partition_length,allele_count,temp_total;
  NumericVector llk(n_partitions);

  arma::umat partition_counts = arma::zeros<arma::umat>(n_snps, 4);
  arma::umat consensus_counts = arma::zeros<arma::umat>(n_snps);

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

//     #calculate log marginal likelihood
//     term1 <- -lgamma(1 + mA+mC+mG+mT)
//     term2 <- lgamma(prior["a", ] + mA) - lgamma(prior["a", ])
//     term2 <- term2 + lgamma(prior["c", ] + mC) - lgamma(prior["c", ])
//     term2 <- term2 + lgamma(prior["g", ] + mG) - lgamma(prior["g", ])
//     term2 <- term2 + lgamma(prior["t", ] + mT) - lgamma(prior["t", ])

  for(p=0; p<n_partitions; p++){
    IntegerVector partition = partitions(p);
    partition_length = partition.size();
    partition_counts.fill(0);

    for (i=0; i<partition_length; i++){
      // fill partition counts matrix and record which snp sites need to be visited.
      arma::sp_mat col(snp_matrix.col(partition(i)-1));
      for (arma::sp_mat::iterator it = col.begin(); it != col.end(); ++it) {
        partition_counts(it.row(), *it-1) += 1;
      }
    }

    arma::sp_umat sparse_partition_counts(partition_counts);

    // calculate llk
    consensus_counts.fill(partition_length);
    llk(p) = term1(partition_length); //-n_snps*pre_lgamma(partition_length+1, 0);

    for (arma::sp_umat::const_iterator it = sparse_partition_counts.begin(); it != sparse_partition_counts.end(); ++it) {
      consensus_counts[it.row()] -= *it;
      // llk(p) += pre_lgamma(*it, prior_index(1+it.col(), it.row()))-pre_lgamma(0, prior_index(1+it.col(), it.row()));
      llk(p) += lgamma(*it + prior(1+it.col(), it.row()))-lgamma(prior(1+it.col(), it.row()));
    }
    for(i=0; i<n_snps; i++){
      // llk(p) += pre_lgamma(consensus_counts[i], prior_index(0, i))-pre_lgamma(0, prior_index(0, i));
      llk(p) += lgamma(consensus_counts[i] + prior(0, i))-lgamma(prior(0, i));
    }

  }

  return(List::create(Named("llk") = llk,
                      Named("term1") = term1));

}
