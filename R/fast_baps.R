#' fast_baps
#'
#' Function to perform Bayesian hierarchical clustering of population structure.
#'
#' [[Rcpp::plugins(cpp11)]]
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param k.init the initial number of clusters to start the bayesian hierarchical clustering from. Defaults to (number of sequences)/4
#' @param n.cores the number of cores to use in clustering
#' @param quiet whether or not to print progress information (default=FALSE)
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' library(Matrix)
#' library(fastbaps)
#' #fasta.file.name <- "../fastbaps_manuscript/data/HIV/hiv_refs_prrt_trim.fas"
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' system.time({baps.hc <- fast_baps(sparse.data)})
#' system.time({baps.hc.2 <- fast_baps(sparse.data, n.cores=4)})
#'
#' @export
fast_baps <- function(sparse.data, k.init=NULL, n.cores=1, quiet=FALSE){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(sparse.data$prior.type %in% c("baps", "mean", "optimised"))) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(n.cores) || n.cores<1) stop("Invalid value for n.cores!")
  if(!is.null(k.init)){
    if(!is.numeric(k.init) | k.init < 0) stop("Invalid value for k.init!")
  }

  n.isolates <- ncol(sparse.data$snp.matrix)
  n.snps <- nrow(sparse.data$snp.matrix)

  if (is.null(k.init)){
    if(sparse.data$prior.type=="baps"){
      k.init <- max(4, ceiling(n.isolates/10))
    } else {
      k.init <- max(4, ceiling(n.isolates/4))
    }
  }

  # Calculate the initial prior paramters
  if(!quiet){
    print("Calculating initial clustering...")
  }
  if(k.init>=n.isolates){
    initial.partition <- as.list(1:n.isolates)
    dk.initial <- rep(0, n.isolates)
  } else {
    h <- get_hclust(sparse.data, quiet)
    if(!quiet){
      print("Calculating initial dk values...")
    }
    dk <- fastbaps:::calc_ddk(sparse.data, h$merge)
    phylo <- fastbaps:::as.phylo.hclust.node.attributes(h, dk$dk[(length(h$labels)+1):length(dk$dk)])
    initial.partition <- split(1:n.isolates, cutree(h, k = min(n.isolates, k.init)))

    dk.initial <- unlist(parallel::mclapply(initial.partition, function(part){
      if(length(part)<=1){
        return(0)
      } else {
        mrca <- ape::getMRCA(phylo, part)
      }
      return(phylo$node.attributes[mrca-length(h$labels)])
    }, mc.cores = n.cores))
  }

  # Cluster the partitions using bhc
  if(!quiet){
    print("Clustering using hierarchical Bayesian clustering...")
  }
  if(n.cores==1){
    baps <- fastbaps:::bhier(sparse.data, initial.partition, dk.initial)
  } else {
    baps <- fastbaps:::bhier_parallel(sparse.data, initial.partition, dk.initial, n.cores)
  }

  hc <- fastbaps:::combine_clusterings(baps, initial.partition)

  hc$labels <- colnames(sparse.data$snp.matrix)


  # Return something similar to the output from hclust along with cuts at different levels.
  return(hc)
}

combine_clusterings <- function(baps, initial.partition){

  merge <- baps$edges
  n.isolates <- length(unlist(initial.partition))
  heights <- rep(0, n.isolates-1)
  group.mems <- -seq_len(n.isolates)
  partition.mems <- rep(0, length(initial.partition))
  new.merge <- matrix(0 , nrow = n.isolates-1, ncol = 2)
  node.count <- 1

  for(i in 1:length(initial.partition)){
    p.len <- length(initial.partition[[i]])
    if(p.len>1){
      new.merge[node.count,1] <- -initial.partition[[i]][[1]]
      new.merge[node.count,2] <- -initial.partition[[i]][[2]]
      heights[[node.count]] <- 0
      node.count <- node.count + 1
      if(p.len>2){
        for(j in 3:(p.len)){
          new.merge[node.count,1] <- node.count - 1
          new.merge[node.count,2] <- -initial.partition[[i]][[j]]
          heights[[node.count]] <- 0
          node.count <- node.count + 1
        }
      }
      partition.mems[[i]] <- node.count - 1
    } else {
      partition.mems[[i]] <- -initial.partition[[i]][[1]]
    }
  }

  row.count <- node.count - 1
  for(i in 1:nrow(merge)){
    if(merge[i,1] < 0){
      new.merge[node.count, 1] <- partition.mems[-merge[i,1]]
    } else {
      new.merge[node.count, 1] <- merge[i,1] + row.count
    }
    if(merge[i,2] < 0){
      new.merge[node.count, 2] <- partition.mems[-merge[i,2]]
    } else {
      new.merge[node.count, 2] <- merge[i,2] + row.count
    }
    heights[[node.count]] <- baps$heights[[i]]
    node.count <- node.count + 1
  }

  hc <- structure(list(merge = new.merge, height = -heights, order = fastbaps:::iorder(new.merge),
                       labels = 1:n.isolates, method = NULL,
                       call = match.call(), dist.method = NULL),
                  class = "hclust")

  return(hc)
}



log_add <- function(u, v){
  if(is.infinite(u) && is.infinite(v)){
    return(-Inf)
  }
  return(max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v))))
}




