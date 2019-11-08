#' root_mll
#'
#' Function to calculate the marginal log likelihood of the root of a hierarchy for use in model comparison.
#'
#' @import Matrix
#'
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param h a hclust object representing the hierarchical clustering
#' @param quiet suppress the printing of extra information (default=FALSE)
#'
#' @return the marginal log likelihood of the root node
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' baps.hc <- fast_baps(sparse.data)
#' root_mll(sparse.data, baps.hc)
#'
#' @export
root_mll <- function(sparse.data, h, quiet=FALSE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(h) %in% c("hclust", "phylo"))) stop("Invalid value for h! Should be a hclust or phylo object!")

  llks <- tree_llk(sparse.data, h$merge)
  root.mll <- llks$ptree[length(llks$ptree)]

  return(root.mll)
}
