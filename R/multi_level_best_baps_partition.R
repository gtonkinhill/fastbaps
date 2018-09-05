#' multi_level_best_baps_partition
#'
#' Function to perform traditional hierarchical clustering, choosing the best partition at each level using the fastbaps approach.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param h a hclust object representing the hierarchical clustering that is to be cut
#' @param levels the number of levels to investigate (default=2)
#' @param n.cores the number of cores to use in clustering
#' @param quiet whether to print additional information or not (default=TRUE)
#'
#' @return a data.frame representing the final clustering at multiple resolutions
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' d <- snp_dist(sparse.data)
#' d <- as.dist(d/max(d))
#' h <- hclust(d, method="ward.D2")
#' multi.res.df <- multi_level_best_baps_partition(sparse.data, h, levels=2)
#'
#' @export
multi_level_best_baps_partition <- function(sparse.data, h, levels=2, n.cores=1, quiet=TRUE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(n.cores) || n.cores<1) stop("Invalid value for n.cores!")
  if(!is.numeric(levels) | levels < 0) stop("Invalid value for levels!")
  if(!(class(h) %in% c("hclust", "phylo"))) stop("Invalid value for h! Should be a hclust or phylo object!")

  if(class(h)=="phylo"){
    if(!ape::is.rooted(h)) stop("phylo object must be rooted")
    h <- ape::multi2di(h)
    nh <- nodeHeights(h)
    tip.edges <- h$edge[,2]<=length(h$tip.label)
    h$edge.length[tip.edges] <- h$edge.length[tip.edges] + (max(nh)-nh[tip.edges,2])
    h <- as.hclust(h)
  }
  if(!all(colnames(sparse.data$snp.matrix) %in% h$labels
  ) || !(all(h$labels %in% colnames(sparse.data$snp.matrix)))){
    stop("Label mismatch between hierarchy and sparse.data!")
  }

  n.isolates <- ncol(sparse.data$snp.matrix)
  n.snps <- nrow(sparse.data$snp.matrix)

  cluster.results <- matrix(1, nrow = n.isolates, ncol = levels+1)

  for (l in seq_len(levels)){
    if(!quiet) print(paste("Clustering at level", l))

    prev.partition <- split(1:n.isolates, cluster.results[,l])

    new.partitions <- rep(NA, n.isolates)
    for (p in seq_along(prev.partition)){
      part <- prev.partition[[p]]
      if (length(part)>4){
        temp.data <- sparse.data
        temp.data$snp.matrix <- temp.data$snp.matrix[,part,drop=FALSE]
        rs <- rowSums(temp.data$snp.matrix>0)
        keep <- rs>0
        keep <- keep & (rs<length(part))
        if(sum(keep)<=1){
          new.partitions[part] <- n.isolates*p*2
          next
        }
        temp.data$snp.matrix <- temp.data$snp.matrix[keep, , drop=FALSE]
        temp.data$prior <- temp.data$prior[, keep, drop=FALSE]
        temp.data$consensus <- temp.data$consensus[keep]

        temp.h <- ape::as.phylo.hclust(h)
        temp.mrca <- ape::getMRCA(temp.h, colnames(temp.data$snp.matrix))
        temp.h <- ape::extract.clade(temp.h, temp.mrca)
        temp.h <- ape::as.hclust.phylo(temp.h)

        new.partitions[part] <- n.isolates*p*2 + fastbaps::best_baps_partition(temp.data, temp.h, quiet = quiet)
      } else {
        new.partitions[part] <- n.isolates*p*2
      }
    }
    #relabel 1:n.clusters
    cluster.results[,l+1] <- as.numeric(factor(new.partitions))

  }

  df <- data.frame(Isolates=colnames(sparse.data$snp.matrix), cluster.results[,2:ncol(cluster.results)],
                   stringsAsFactors = FALSE)
  colnames(df)[2:ncol(df)] <- paste("Level", 1:levels)

  return(df)
}
