#' optimise_prior
#'
#' Function to optimise the prior using Bayes Factors and grid search
#'
#' @import Matrix
#'
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param grid.interval the upper and lower bound for the hyperparameter to optimise over (default=c(5e-4, 10))
#' @param type one of "flat" or "hc" indicating a flat prior or the approach of Heller et al.
#' @param n.cores number of cores to use (currently not implemented)
#'
#' @return a sparse.data object with the prior optimised via grid search
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' sparse.data <- optimise_prior(sparse.data)
#'
#' @export
optimise_prior <- function(sparse.data, grid.interval=c(5e-4, 10), type = "symmetric", n.cores=1){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(type %in% c("symmetric", "hc", "optimise.baps", "baps", "ref"))) stop("Invalid value for type. Must be one of 'symmetric', 'hc', 'ref' or 'baps'")
  if(!all(grid.interval>0)) stop("grid values must greater than 0")

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
    initial.prior <- ceiling(initial.prior*1000)/1000
  } else if (type=="optimise.baps") {
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), nrow(sparse.data$snp.matrix)),
                              rowSums(sparse.data$snp.matrix==1),
                              rowSums(sparse.data$snp.matrix==2),
                              rowSums(sparse.data$snp.matrix==3),
                              rowSums(sparse.data$snp.matrix==4)), nrow = 5, byrow = TRUE)
    initial.prior[1,] <- initial.prior[1,] - colSums(initial.prior[2:5,])
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior <- initial.prior>0
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
  } else if (type=="baps") {
    initial.prior <- matrix(c(rep(ncol(sparse.data$snp.matrix), nrow(sparse.data$snp.matrix)),
                              rowSums(sparse.data$snp.matrix==1),
                              rowSums(sparse.data$snp.matrix==2),
                              rowSums(sparse.data$snp.matrix==3),
                              rowSums(sparse.data$snp.matrix==4)), nrow = 5, byrow = TRUE)
    initial.prior[1,] <- initial.prior[1,] - colSums(initial.prior[2:5,])
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
    initial.prior <- initial.prior>0
    initial.prior <- t(t(initial.prior)/colSums(initial.prior))
  } else if (type=="ref") {
    initial.prior <- matrix(1, nrow = nrow(sparse.data$prior), ncol = ncol(sparse.data$prior))
    initial.prior[1,] <- 2
  } else {
    initial.prior <- matrix(1, nrow = nrow(sparse.data$prior), ncol = ncol(sparse.data$prior))
  }
  sparse.data$prior <- initial.prior

  if (type=="baps") return(sparse.data)

  opt <- optimise(calc_prior_prob, grid.interval, sparse.data, initial.prior, h, maximum = TRUE, tol=1e-3)
  cc <- round(opt$maximum, digits = 3)

  sparse.data$prior <- initial.prior*cc
  sparse.data$prior[sparse.data$prior<1e-3] <- 1e-3

  ##Check hyperparameter is not too close to grid edges
  if (any(abs(grid.interval-cc)<5e-3)){
    warning("Inferred hyperparameter is very close to interval boundries! Consider changing the interval.")
  }

  print(paste("Optimised hyperparameter:", cc))

  sparse.data$prior.type = "optimised"

  return(sparse.data)
}

calc_prior_prob <- function(cc, temp.sparse.data, temp.initial.prior, temp.h){
  cc <- round(cc, digits = 3) #as otherwise we run into issue with rounding in the llk calculation due to the way lgamma is calculated
  temp.sparse.data$prior <- temp.initial.prior * cc
  temp.sparse.data$prior[temp.sparse.data$prior<1e-3] <- 1e-3
  llks <- fastbaps:::tree_llk(temp.sparse.data, temp.h$merge)
  prob.tree <- llks$ptree[length(llks$ptree)]
  return(prob.tree)
}


