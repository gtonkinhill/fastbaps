#' calc_marginal_llk
#'
#' Function to calculate the log marginal likelihood of a given partition
#'
#' @import Matrix
#'
#' @param sparse.matrix a sparse SNP data object returned from import_fasta_sparse_nt
#' @param partition the partition for which the log marginal likelihood is to be calculated
#'
#' @return the log marginal likelihood
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' partition <- 1:ncol(sparse.data$snp.matrix)
#' calc_marginal_llk(sparse.data, partition)
#'
#' @export
calc_marginal_llk <- function(sparse.data, partition){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(partition)) stop("Invalid value for partition! Should be a numeric vector with length equal to number of sequences!")
  if(ncol(sparse.data$snp.matrix)!=length(partition)) stop("Invalid value for partition! Should be a numeric vector with length equal to number of sequences!")
  if(!all(c(1:length(unique(partition))) %in% partition)) stop("partition should be a numeric vector consisting of the integers 1 to n.clusters!")

  partition.list <- split(1:ncol(sparse.data$snp.matrix), partition)

  mllk <- fastbaps:::part_llks(sparse.data, partition.list)$llk

  return(sum(mllk))
}
