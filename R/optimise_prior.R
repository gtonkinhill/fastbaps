#' optimise_prior
#'
#' Function to optimise the prior using Bayes Factors and grid search
#'
#' @import Matrix
#'
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param grid.vals the grid values of the hyperparameter to optimise over (default=c(0:20/5))
#' @param type one of "flat" or "hc" indicating a flat prior or the approach of Heller et al.
#' @param n.cores number of cores to use
#'
#' @return a sparse.data object with the prior optimised via grid search
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' sparse.data <- optimise_prior(sparse.data)
#'
#' @export
optimise_prior <- function(sparse.data, grid.vals=c(1:100/20), type = "flat", n.cores=1){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(type %in% c("flat", "hc", "combined"))) stop("Invalid value for type. Must be one of 'flat' or 'hc'")
  if(!all(grid.vals>0)) stop("grid values must greater than 0")

  MIN_RES <- 1e-3

  #create hclust object
  h <- get_hclust(sparse.data, TRUE)

  if(type=="hc"){
    #initialise prior
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), nrow(sparse.data$snp.matrix)),
                              rowSums(sparse.data$snp.matrix==1),
                              rowSums(sparse.data$snp.matrix==2),
                              rowSums(sparse.data$snp.matrix==3),
                              rowSums(sparse.data$snp.matrix==4)), nrow = 5, byrow = TRUE)
    initial.prior[1,] <- initial.prior[1,] - colSums(initial.prior[2:5,])
    initial.prior <- initial.prior + 1
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior[initial.prior<MIN_RES] <- MIN_RES
    sparse.data$prior <- initial.prior

    ptrees <- unlist(parallel::mclapply(grid.vals, function(cc){
      sparse.data$prior <- initial.prior * cc
      sparse.data$prior[sparse.data$prior<MIN_RES] <- MIN_RES
      llks <- fastbaps:::tree_llk(sparse.data, h$merge)
      return(llks$ptree[length(llks$ptree)])
    }, mc.cores = n.cores))

    if(!all(ptrees<0)) stop("Error: tree probabilities exceed 1!")

    bfs <- ptrees[[1]] - ptrees

    cc <- grid.vals[which.min(bfs)]
    sparse.data$prior <- initial.prior*cc
  } else if (type=="combined"){
    #initialise prior
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), nrow(sparse.data$snp.matrix)),
                              rowSums(sparse.data$snp.matrix==1),
                              rowSums(sparse.data$snp.matrix==2),
                              rowSums(sparse.data$snp.matrix==3),
                              rowSums(sparse.data$snp.matrix==4)), nrow = 5, byrow = TRUE)
    initial.prior[1,] <- initial.prior[1,] - colSums(initial.prior[2:5,])
    initial.prior <- initial.prior + 1
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior[initial.prior<MIN_RES] <- MIN_RES
    sparse.data$prior <- initial.prior

    ptrees <- unlist(parallel::mclapply(grid.vals, function(cc){
      sparse.data$prior <- initial.prior + cc
      sparse.data$prior[sparse.data$prior<MIN_RES] <- MIN_RES
      llks <- fastbaps:::tree_llk(sparse.data, h$merge)
      return(llks$ptree[length(llks$ptree)])
    }, mc.cores = n.cores))

    bfs <- ptrees[[1]] - ptrees

    cc <- grid.vals[which.min(bfs)]
    sparse.data$prior <- initial.prior + cc
  } else {
    initial.prior <- matrix(MIN_RES, nrow = nrow(sparse.data$prior), ncol = ncol(sparse.data$prior))
    sparse.data$prior <- initial.prior

    ptrees <- unlist(parallel::mclapply(grid.vals, function(cc){
      sparse.data$prior <- initial.prior + cc
      sparse.data$prior[sparse.data$prior<MIN_RES] <- MIN_RES
      sparse.data$prior <- t(t(sparse.data$prior)/colSums(sparse.data$prior))
      llks <- fastbaps:::tree_llk(sparse.data, h$merge)
      return(llks$ptree[length(llks$ptree)])
    }, mc.cores = n.cores))

    bfs <- ptrees[[1]] - ptrees

    # mllks <- fastbaps:::compare_prior_grid(sparse.data, grid.vals)
    # bfs <- mllks[[1]] - mllks

    cc <- grid.vals[which.min(bfs)]
    sparse.data$prior <- initial.prior+cc
  }

  print(paste("Optimised hyperparameter:", cc))

  sparse.data$prior.type = "optimised"

  return(sparse.data)
}
