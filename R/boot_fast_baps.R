#' boot_fast_baps
#'
#' Function to perform bootstrap replicated of fastbaps.
#'
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param k.init the initial number of clusters to start the bayesian hierarchical clustering from. Defaults to (number of sequences)/4
#' @param n.replicates the number of bootstrap replicates to perform (default=100)
#' @param hc.method the type of initial hierarchical clustering to use. Can be with 'ward' or 'genie' (default='ward')
#' @param n.cores the number of cores to use in clustering
#' @param quiet whether or not to print progress information (default=FALSE)
#'
#' @return a co-occurence count matrix
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' boot.result <- boot_fast_baps(sparse.data)
#' dendro <- as.dendrogram(fast_baps(sparse.data))
#' heatmap(boot.result, dendro, dendro)
#'
#'
#' @export
boot_fast_baps <- function(sparse.data, k.init=NULL, n.replicates=100, hc.method='ward',
                           n.cores=1, quiet=TRUE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(sparse.data$prior.type %in% c("baps", "mean", "optimised"))) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(n.cores) || n.cores<1) stop("Invalid value for n.cores!")
  if(!is.null(k.init)){
    if(!is.numeric(k.init) | k.init < 0) stop("Invalid value for k.init!")
  }
  if(!is.numeric(n.replicates) | n.replicates < 1) stop("Invalid value for replicates!")
  if(!(hc.method %in% c("ward", "genie"))) stop("Invalid hc.method!")

  n.isolates <- ncol(sparse.data$snp.matrix)
  n.snps <- nrow(sparse.data$snp.matrix)

  boot.matrix <- matrix(0, nrow = n.isolates, ncol = n.isolates)
  colnames(boot.matrix) <- colnames(sparse.data$snp.matrix)
  rownames(boot.matrix) <- colnames(sparse.data$snp.matrix)

  for (r in 1:n.replicates){
    boot.sample <- sample(1:n.snps, n.snps, replace = TRUE)
    temp.data <- sparse.data
    temp.data$snp.matrix <- temp.data$snp.matrix[boot.sample,]
    temp.data$prior <- temp.data$prior[,boot.sample]
    temp.data$consensus <- temp.data$consensus[boot.sample]
    temp.hc <- fast_baps(temp.data, k.init, hc.method, n.cores, quiet)
    temp.clusters <- best_baps_partition(temp.data, temp.hc, quiet)
    for (c in 1:max(temp.clusters)){
      if(sum(temp.clusters==c)<=1) next
      pairs <- t(utils::combn(names(temp.clusters)[temp.clusters==c], 2))
      boot.matrix[pairs] = boot.matrix[pairs, drop=FALSE] + 1
      boot.matrix[pairs[,c(2,1), drop=FALSE]] = boot.matrix[pairs[,c(2,1), drop=FALSE], drop=FALSE] + 1
    }
  }
  diag(boot.matrix) <- n.replicates

  # Return a matrix of co-occurence
  return(boot.matrix)
}
