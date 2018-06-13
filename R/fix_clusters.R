#' fix_clusters
#'
#' Function to iteratively update an intially clustering, greedily maximising the BAPS marginal likelihood.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param clusters the initial clustering to be improved. Given as a vector of cluster memberships of length equal to the number of sequences
#' @param n.iterations the number of greedy maximisation iterations to perform
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' sim.matrix <- snp_similarity(sparse.data)
#' x <- as.dist(1-sim.matrix/max(sim.matrix))
#' phylo <- ape::as.phylo(hclust(x, method="average"))
#' clusters <- best_baps_partition(sparse.data, phylo)
#' clusters <- fix_clusters(sparse.data, clusters)
#'
#' @export
fix_clusters <- function(sparse.data, clusters, n.iterations=1){

  n.isolates <- ncol(sparse.data$snp.matrix)
  n.snps <- nrow(sparse.data$snp.matrix)
  partitions <- split(1:n.isolates, clusters)
  current.llks <- fastbaps:::part_llks(sparse.data, partitions)$llk

  is.improved <- FALSE
  for (it in 1:n.iterations){
    for (i in 1:n.isolates){
      #find best potential swap
      temp.partitions <- lapply(partitions, function(p) unique(unlist(c(i, p))))
      temp.part.llks <- fastbaps:::part_llks(sparse.data, temp.partitions)$llk
      diff <- current.llks-temp.part.llks
      diff[[clusters[i]]] <- Inf
      best.move <- which.min(diff)

      #recalibrate temp partition
      temp.clusters <- clusters
      temp.clusters[[i]] <- best.move
      temp.clusters <- as.numeric(factor(temp.clusters)) #re-number in case a cluster has been removed.
      temp.partitions <- split(1:n.isolates, temp.clusters)
      temp.part.llks <- fastbaps:::part_llks(sparse.data, temp.partitions)$llk

      #test if it improves things
      if(sum(current.llks) < sum(temp.part.llks)){
        current.llks <- temp.part.llks
        partitions <- temp.partitions
        clusters <- temp.clusters
        is.improved <- TRUE
        print("Improved!")
      }
    }
  }

  return(clusters)
}
