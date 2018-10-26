#' import_fasta_sparse_nt
#'
#' Imports a fasta file to a sparse matrix representing SNPs from the consensus
#'
#' @useDynLib fastbaps
#' @importFrom Rcpp sourceCpp
#' @import Matrix
#'
#' @param fasta path to the fasta file or an ape DNAbin object
#' @param prior the type of prior to use. Can be one of 'baps' and 'mean' (default=baps)
#' @param check.fasta whether to check the fasta file for issue. Slows things down a little but may be avoided for very large fasta files.
#'
#' @return A sparse matrix reprsentation of the SNPs (different to the consensus sequence)
#'
#' @examples
#' fasta <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' fasta <- ape::read.FASTA(fasta)
#' sparse.data <- import_fasta_sparse_nt(fasta)
#'
#' @export
import_fasta_sparse_nt <- function(fasta, prior='baps', check.fasta=TRUE){

  # Check inputs
  if(class(fasta)!="DNAbin") {
    if(!file.exists(fasta)) stop(paste("Can't locate file", fasta))
    if(!is.logical(check.fasta)) stop("check.fasta should be one of TRUE/FALSE!")
    if(check.fasta){
      # Cheat a bit by checking the file using ape
      invisible(utils::capture.output(ape::read.FASTA(fasta, type = "DNA")))
    }
  }
  if(!is.character(prior)) stop("Invalid input for prior parameter!")
  if(!(prior %in% c("baps", "mean"))) stop("Invalid input for prior parameter!")

  if(class(fasta)=="DNAbin"){
    fasta <- as.character(as.matrix(fasta))
    ij <- which(t(fasta) != fasta[1,], arr.ind = TRUE)
    fasta[fasta=='a'] <- 1
    fasta[fasta=='c'] <- 2
    fasta[fasta=='g'] <- 3
    fasta[fasta=='t'] <- 4
    fasta[fasta=='-'] <- 5
    fasta[fasta=='n'] <- 5
    fasta <- apply(fasta, 2, as.numeric)
    snp.data <- list(num.seqs=nrow(fasta),
                     consensus=fasta[1,]-1,
                     seq.length=ncol(fasta),
                     seq.names=rownames(fasta))

    snp.matrix <- t(sparseMatrix(i=ij[,1], j=ij[,2], x=t(fasta)[ij],
                                 dims = c(snp.data$seq.length, snp.data$num.seqs),
                                 dimnames = list(1:snp.data$seq.length, snp.data$seq.names)))

  } else {
    snp.data <- import_fasta_to_vector_each_nt(fasta)
    snp.data$seq.names <-  gsub("^>","",snp.data$seq.names)

    snp.matrix <- sparseMatrix(i=snp.data$i,
                               j=snp.data$j,
                               x=snp.data$x,
                               dims = c(snp.data$num.seqs, snp.data$seq.length),
                               dimnames = list(snp.data$seq.names, 1:snp.data$seq.length))
  }

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
