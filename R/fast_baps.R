#' fast_baps
#'
#' Function to perform Bayesian hierarchical clustering of population structure.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param k.init the initial number of clusters to start the bayesian hierarchical clustering from. Defaults to (number of sequences)/4
#' @param n.cores the number of cores to use in clustering
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' system.time({baps.hc <- fast_baps(sparse.data)})
#' system.time({baps.hc.2 <- fast_baps(sparse.data, n.cores=8)})
#' baps.hc$height-baps.hc.2$height
#'
#' @export
fast_baps <- function(sparse.data, k.init=NULL, n.cores=1){

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
  if(k.init>=n.isolates){
    initial.partition <- as.list(1:n.isolates)
    dk.initial <- rep(0, n.isolates)
  } else {
    snp.dist <- as.matrix(tcrossprod(t(sparse.data$snp.matrix>0)))
    snp.dist <- as.dist((max(snp.dist)-snp.dist)/max(snp.dist))
    h <- hclust(snp.dist, method = "ward.D2")
    phylo <- ape::as.phylo(h)
    initial.partition <- split(1:n.isolates, cutree(h, k = min(n.isolates, k.init)))
    dk.initial <- unlist(lapply(initial.partition, function(part){
      if(length(part)<=1){
        return(0)
      }
      mrca <- ape::getMRCA(phylo, part)
      temp.clade <- ape::extract.clade(phylo, mrca)
      return(fastbaps:::get_dk(temp.clade))
    }))
    dk.initial <- rep(0, length(initial.partition))
  }

  # Cluster the partitions using bhc
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

get_dk <- function(clade){
  ntip <- length(clade$tip.label)
  if(ntip<=2){
    return(0)
  } else {
    children <- phangorn::Children(clade, clade$edge[1,1])
    if (children[[1]] <= ntip){
      dl <- 0
    } else {
      dl <- get_dk(ape::extract.clade(clade, children[[1]]))
    }
    if (children[[2]] <= ntip){
      dr <- 0
    } else {
      dr <- get_dk(ape::extract.clade(clade, children[[2]]))
    }
    return(log_add(lgamma(ntip), dl+dr))
  }
}

log_add <- function(u, v){
  if(is.infinite(u) && is.infinite(v)){
    return(-Inf)
  }
  return(max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v))))
}



