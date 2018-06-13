#' snp_similarity
#'
#' Function to calculate pairwise snp similarity matrix
#'
#' @useDynLib fastbaps
#' @importFrom Rcpp sourceCpp
#' @import Matrix
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

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")

  shared.snps <- as.matrix(tcrossprod(t(sparse.data$snp.matrix==1)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==2)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==3)))
  shared.snps <- shared.snps + as.matrix(tcrossprod(t(sparse.data$snp.matrix==4)))

  return(shared.snps)
}
