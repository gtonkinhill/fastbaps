#' as.phylo.hclust.node.attributes
#' get_hclust
#'
#' Function to generate hclust object. Uses initial kmeans if dataset is too large.
#'
#' @import ape
#' @import Matrix
#'
#' @param sparse.data sparse data object generated using
#' @param quiet whether to supress additional printing information
#' @param method which hierarchical clustering method. One of 'ward' or 'genie'. (default='ward')
#' @param n.cores n.cores to use (not implemented)
#'
#'
# adpated from ape as.phylo
get_hclust <- function(sparse.data, quiet, method="ward", n.cores=1){

  MAX_CLUSTER_TO_GENIE <- 10000
  n.isolates <- ncol(sparse.data$snp.matrix)

  stopifnot(method %in% c("ward", "genie"))


  if(n.isolates<MAX_CLUSTER_TO_GENIE){
    snp.dist <- as.matrix(tcrossprod(t(sparse.data$snp.matrix>0)))
    snp.dist <- stats::as.dist((max(snp.dist)-snp.dist)/max(snp.dist))
    if(method=="genie"){
      h <- genie::hclust2(d=snp.dist, useVpTree=FALSE)
      h$labels <- colnames(sparse.data$snp.matrix)
    } else {
      h <- stats::hclust(snp.dist, method = "ward.D2")
    }
  } else {
    if(!quiet){
      print("Large number of sequences so using an initial PCA and the genie hierarchical clustering algorithm.")
    }
    temp.matrix <- 1*t(sparse.data$snp.matrix>0)
    # pc <- irlba::prcomp_irlba(temp.matrix, n=50)
    svd <- irlba::irlba(temp.matrix, nv=50, tol=0.1, center=colMeans(temp.matrix))
    x <- t(t(svd$u) * svd$d)
    if(method=="genie"){
      h <- genie::hclust2(d="euclidean", objects = x, useVpTree=FALSE, thresholdGini = 1)
      h$labels <- colnames(sparse.data$snp.matrix)
    } else {
      h <- fastcluster::hclust.vector(X = x, method = "ward")
      h$labels <- colnames(sparse.data$snp.matrix)
    }
    gc()

  }
  return(h)

}
