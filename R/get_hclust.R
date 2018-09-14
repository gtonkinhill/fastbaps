#' as.phylo.hclust.node.attributes
#' get_hclust
#'
#' Function to generate hclust object. Uses initial kmeans if dataset is too large.
#'
#' @import ape
#' @import Matrix
#'
#' @param sparse.data
#' @param attribute
#' @param n.cores
#'
#'
# adpated from ape as.phylo
get_hclust <- function(sparse.data, quiet, n.cores=1){

  MAX_CLUSTER_TO_GENIE <- 10000
  n.isolates <- ncol(sparse.data$snp.matrix)

  if(n.isolates<MAX_CLUSTER_TO_GENIE){
    snp.dist <- as.matrix(tcrossprod(t(sparse.data$snp.matrix>0)))
    snp.dist <- as.dist((max(snp.dist)-snp.dist)/max(snp.dist))
    h <- stats::hclust(snp.dist, method = "ward.D2")
  } else {
    if(!quiet){
      print("Large number of sequences so using an initial PCA and the genie hierarchical clustering algorithm.")
    }
    temp.matrix <- 1*t(sparse.data$snp.matrix>0)
    svd <- irlba::irlba(temp.matrix, nv=50, tol=0.1, center=colMeans(temp.matrix))
    h <- genie::hclust2(d="euclidean", objects = svd$u, useVpTree=FALSE, thresholdGini = 1)
    h$labels <- colnames(sparse.data$snp.matrix)
    gc()

  }
  return(h)

}
