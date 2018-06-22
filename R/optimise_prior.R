#' optimise_prior
#'
#' Function to optimise the prior using Bayes Factors and grid search
#'
#' @import Matrix
#'
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param grid.vals the grid values of the hyperparameter to optimise over (default=c(0:20/5))
#'
#' @return a sparse.data object with the prior optimised via grid search
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' sparse.data <- optimise_prior(sparse.data)
#'
#' @export
optimise_prior <- function(sparse.data, grid.vals=c(0:1000/200)){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")

  n.isolates <- ncol(sparse.data$snp.matrix)

  partition <- list(1:n.isolates)

  initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), nrow(sparse.data$snp.matrix)),
                           rowSums(sparse.data$snp.matrix==1),
                           rowSums(sparse.data$snp.matrix==2),
                           rowSums(sparse.data$snp.matrix==3),
                           rowSums(sparse.data$snp.matrix==4)), nrow = 5, byrow = TRUE)
  initial.prior[1,] <- initial.prior[1,] - colSums(initial.prior[2:5,])
  initial.prior <- t(t(initial.prior)/colSums(initial.prior))
  sparse.data$prior <- initial.prior

  mllks <- fastbaps:::compare_prior_grid(sparse.data, grid.vals)
  bfs <- mllks[[1]] - mllks

  cc <- grid.vals[which.min(bfs)]
  sparse.data$prior <- initial.prior+cc
  sparse.data$prior <- t(t(sparse.data$prior)/colSums(sparse.data$prior))

  # grid.vals <- 1:200/10
  # bfs <- unlist(lapply(grid.vals, function(cc){
  #   sparse.data$prior <- initial.prior*cc
  #   mllk2 <- fastbaps:::part_llks(sparse.data, partition)$llk
  #   return(mllk1-mllk2)
  # }))
  # plot(bfs)
  #
  # cc <- grid.vals[which.min(bfs)]
  # sparse.data$prior <- initial.prior*cc

  print(paste("Optimised hyperparameter:", cc))

  sparse.data$prior.type = "optimised"

  return(sparse.data)
}
