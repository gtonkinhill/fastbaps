#' snp_similarity
#'
#' Function to calculate pairwise snp similarity matrix
#'
#' @useDynLib fastbaps
#' @importFrom Rcpp sourceCpp
#'
#' @param sparse.matrix a sparse SNP data object returned from import_fasta_sparse_nt
#'
#' @return A pairwise snp similarity matrix
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' snp_similarity(sparse.data)
#'
#' @export
snp_similarity <- function(sparse.data){

  shared.snps <- as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==1)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==2)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==3)))
  shared.snps <- shared.snps + as.matrix(Matrix::tcrossprod(Matrix::t(sparse.data$snp.matrix==4)))

  return(shared.snps)
}
