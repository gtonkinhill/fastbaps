#' snp_dist
#'
#' Function to calculate pairwise snp distance matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#'
#' @return A pairwise snp distance matrix
#'
#' @examples
#' fasta.file <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' snp_dist(sparse.data)
#'
#' @export
snp_dist <- function(sparse.data){

  n.isolates <- ncol(sparse.data$snp.matrix)

  shared.snps <- as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==1)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==2)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==3)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==4)))

  total.snps <- Matrix::colSums(sparse.data$snp.matrix>0)
  differing.snps <- (matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = TRUE) +
                       matrix(rep(total.snps, n.isolates), nrow = n.isolates, byrow = FALSE) -
                       as.matrix(tcrossprod(t(snp.data$snp.matrix>0))) -
                       shared.snps)

  return(differing.snps)
}
