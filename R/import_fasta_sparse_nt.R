#' import_fasta_sparse_nt
#'
#' Imports a fasta file to a sparse matrix representing SNPs from the consensus
#'
#' @useDynLib fastbaps, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#'
#' @param fasta.file.name path to the fasta file
#'
#' @return A sparse matrix reprsentation of the SNPs (different to the consensus sequence)
#'
#' @examples
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#'
#' @export
import_fasta_sparse_nt <- function(fasta.file.name){

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
  prior <- matrix(c(rep(nrow(snp.matrix), ncol(snp.matrix)),
                    colSums(snp.matrix==1),
                    colSums(snp.matrix==2),
                    colSums(snp.matrix==3),
                    colSums(snp.matrix==4)), nrow = 5, byrow = TRUE)
  prior[1,] <- prior[1,]- colSums(prior[2:5,])
  prior <- t(t(prior)/colSums(prior))


  return(list(snp.matrix=t(snp.matrix), consensus=snp.data$consensus,
              prior=prior))
}
