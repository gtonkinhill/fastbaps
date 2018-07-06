#' as.phylo.hclust.node.attributes
#' get_hclust
#'
#' Function to generate hclust object. Uses initial kmeans if dataset is too large.
#'
#' @import ape
#'
#' @param sparse.data
#' @param attribute
#'
#'
# adpated from ape as.phylo
get_hclust <- function(sparse.data, quiet){

  MAX_CLUSTER_TO_KMEANS <- 10000
  n.isolates <- ncol(sparse.data$snp.matrix)

  if(n.isolates<MAX_CLUSTER_TO_KMEANS){
    snp.dist <- as.matrix(tcrossprod(t(sparse.data$snp.matrix>0)))
    snp.dist <- as.dist((max(snp.dist)-snp.dist)/max(snp.dist))
    h <- stats::hclust(snp.dist, method = "ward.D2")
  } else {
    if(!quiet){
      print("Large number of sequences so doing an initial split using kmeans...")
    }
    #as we are only interested in the lowere branches of the dendrogram we can take advatange of a fast
    #algorithm like kmeans to reduce the complexity of the problem
    pc <- irlba::prcomp_irlba(1*t(sparse.data$snp.matrix>0), n=50)
    k.n.clusters <- ceiling(n.isolates/MAX_CLUSTER_TO_KMEANS)
    k <- kmeans(pc$x, centers = k.n.clusters, iter.max = 50, nstart = 5)$cluster
    while(max(table(k))>MAX_CLUSTER_TO_KMEANS){
      k.n.clusters <- k.n.clusters*2
      k <- kmeans(pc$x, centers = k.n.clusters, iter.max = 50, nstart = 5)$cluster
    }

    hclist <- lapply(split(1:n.isolates, k), function(k.part){
      snp.dist <- as.matrix(tcrossprod(t(sparse.data$snp.matrix[,k.part]>0)))
      snp.dist <- as.dist((max(snp.dist)-snp.dist)/max(snp.dist))
      ht <- hclust(snp.dist, method = "ward.D2")
      ht$height <- round(ht$height, 6) #to stop us running into precision issues
      return(ht)
    })
    h <- merge_hclust(hclist)
  }
  return(h)

}


merge_hclust <- function(hclist) {
  d <- as.dendrogram(hclist[[1]])
  for (i in 2:length(hclist)) {
    # print(i)
    d <- merge(d, as.dendrogram(hclist[[i]]))
  }
  as.hclust(d)
}
