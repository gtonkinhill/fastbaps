#' best_baps_partition
#'
#' Function to combine smaller clusters from a fast hierarchical algorithm to maximise the BAPS likelihood.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param h a hclust object representing the hierarchical clustering that is to be cut
#' @param quiet
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' d <- snp_dist(sparse.data)
#' d <- as.dist(d/max(d))
#' h <- hclust(d, method="ward.D2")
#' partition <- best_baps_partition(sparse.data, h)
#'
#' newick.file.name <- system.file("extdata", "seqs.fa.treefile", package = "fastbaps")
#' iqtree <- phytools::read.newick(newick.file.name)
#' h <- phytools::midpoint.root(iqtree)
#' best.partition <- best_baps_partition(sparse.data, h)
#'
#' @export
best_baps_partition <- function(sparse.data, h, quiet=FALSE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(h) %in% c("hclust", "phylo"))) stop("Invalid value for h! Should be a hclust or phylo object!")

  if(class(h)=="phylo"){
    if(!ape::is.rooted(h)) stop("phylo object must be rooted")
    h <- as.hclust(phytools::force.ultrametric(h, method = "extend"))
  }
  if(!all(colnames(sparse.data$snp.matrix) %in% h$labels
          ) || !(all(h$labels %in% colnames(sparse.data$snp.matrix)))){
    stop("Label mismatch between hierarchy and sparse.data!")
  }

  sparse.data$snp.matrix <- sparse.data$snp.matrix[,match(h$labels, colnames(sparse.data$snp.matrix))]

  if(!quiet){
    print("Calculating node marginal llks...")
  }
  llks <- fastbaps:::tree_llk(sparse.data, h$merge)
  n.isolates <- ncol(sparse.data$snp.matrix)

  if(!quiet){
    print("Finding best partition...")
  }

  threshold <- log(0.5)
  rk <- llks$rk[(n.isolates+1):length(llks$rk)]
  clusters <- -c(1:n.isolates)
  for (i in 1:nrow(h$merge)){
    if(rk[[i]]>threshold){
      clusters[clusters %in% h$merge[i,]] <- i
    }
  }
  clusters <- as.numeric(factor(clusters))
  names(clusters) <- h$labels

  return(clusters)
}
