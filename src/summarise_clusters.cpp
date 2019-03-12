// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
IntegerVector summarise_clusters(NumericMatrix merge, NumericVector rk, double threshold, int n_isolates) {

  int i;
  NumericVector temp_rk;
  temp_rk = rk[Rcpp::Range(0, (merge.nrow()-1))];

  for (i=(merge.nrow()-1); i>=0; i--){
    if(temp_rk(i)>threshold){
      if (merge(i,0)>0){
        temp_rk(merge(i,0)-1) = threshold+1;
      }
      if (merge(i,1)>0){
        temp_rk(merge(i,1)-1) = threshold+1;
      }
    }
  }


  std::vector< std::vector<int> > clusters(merge.nrow());
  IntegerVector keep(merge.nrow(), 0);


  int j;

  for (i=0; i<merge.nrow(); i++){
    if(temp_rk(i)>threshold){
      if(merge(i,0)<0){
        clusters[i].push_back(abs(merge(i,0)));
      } else {
        clusters[i].insert(clusters[i].end(), clusters[merge(i,0)-1].begin(), clusters[merge(i,0)-1].end());
        keep(merge(i,0)-1) = 1;
        clusters[merge(i,0)-1].clear();
      }
      if(merge(i,1)<0){
        clusters[i].push_back(abs(merge(i,1)));
      } else {
        clusters[i].insert(clusters[i].end(), clusters[merge(i,1)-1].begin(), clusters[merge(i,1)-1].end());
        keep(merge(i,1)-1) = 1;
        clusters[merge(i,1)-1].clear();
      }
    } else {
      keep(i) = 1;
    }
  }

  IntegerVector final_clust(n_isolates, 0);
  int n_clust=0;
  for (i=0; i<clusters.size(); i++){
    n_clust+=1;
    if (keep(i)==0){
      for (j=0; j<clusters[i].size(); j++){
        final_clust(clusters[i][j]-1) = i+1;
      }
    }
  }

  for (i=0; i< n_isolates; i++){
    if(final_clust(i)==0){
      n_clust+=1;
      final_clust(i)=n_clust;
    }
  }

  return(final_clust);

}
