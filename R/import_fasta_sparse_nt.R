#' import_fasta_sparse_nt
#'
#' Imports a fasta file to a sparse matrix representing SNPs from the consensus
#'
#' @useDynLib fastbaps
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#'
#' @param fasta.file.name path to the fasta file
#' @param prior the type of prior to use. Can be one of 'baps' and 'mean' (default=baps)
#'
#' @return A sparse matrix reprsentation of the SNPs (different to the consensus sequence)
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#'
#' @export
import_fasta_sparse_nt <- function(fasta.file.name, prior='baps'){

  # Check inputs
  if(!file.exists(fasta.file.name)) stop(paste("Can't locate file", fasta.file.name))
  # Cheat a bit by checking the file using ape
  invisible(capture.output(ape::read.FASTA(fasta.file.name, type = "DNA")))
  if(!is.character(prior)) stop("Invalid input for prior parameter!")
  if(!(prior %in% c("baps", "mean"))) stop("Invalid input for prior parameter!")

  snp.data <- fastbaps:::import_fasta_to_vector_each_nt(fasta.file.name)
  snp.data$seq.names <-  gsub("^>","",snp.data$seq.names)

  snp.matrix <- sparseMatrix(i=snp.data$i,
                             j=snp.data$j,
                             x=snp.data$x,
                             dims = c(snp.data$num.seqs, snp.data$seq.length),
                             dimnames = list(snp.data$seq.names, 1:snp.data$seq.length))

  #Remove columns where the consensus is gapped
  snp.matrix <- snp.matrix[, snp.data$consensus!=4]
  snp.data$consensus <- snp.data$consensus[snp.data$consensus!=4]


  #remove conserved columns
  conserved <- colSums(snp.matrix>0) == 0
  snp.matrix <- snp.matrix[,!conserved]
  snp.data$consensus <- snp.data$consensus[!conserved]


  #prior matrix (number of unique allels at each locus/denominator)
  prior.matrix <- matrix(c(rep(nrow(snp.matrix), ncol(snp.matrix)),
                    colSums(snp.matrix==1),
                    colSums(snp.matrix==2),
                    colSums(snp.matrix==3),
                    colSums(snp.matrix==4)), nrow = 5, byrow = TRUE)
  prior.matrix[1,] <- prior.matrix[1,]- colSums(prior.matrix[2:5,])
  prior.matrix <- t(t(prior.matrix)/colSums(prior.matrix))

  if (prior=="baps"){
    prior.matrix <- prior.matrix>0
    prior.matrix <- t(t(prior.matrix)/colSums(prior.matrix))
  }


  return(list(snp.matrix=t(snp.matrix), consensus=snp.data$consensus,
              prior=prior.matrix, prior.type=prior))
}
