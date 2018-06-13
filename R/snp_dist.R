#' snp_dist
#'
#' Function to calculate pairwise snp distance matrix
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#'
#' @return A pairwise snp distance matrix
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' snp_dist(sparse.data)
#'
#' @export
snp_dist <- function(sparse.data){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")

  n.isolates <- ncol(sparse.data$snp.matrix)

  shared.snps <- as.matrix(tcrossprod(t(sparse.data$snp.matrix==1)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==2)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==3)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==4)))

  total.snps <- colSums(sparse.data$snp.matrix>0)
  differing.snps <- (matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = TRUE) +
                       matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = FALSE) -
                       as.matrix(tcrossprod(t(sparse.data$snp.matrix>0))) -
                       shared.snps)

  return(differing.snps)
}
