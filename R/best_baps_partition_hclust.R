#' best_baps_partition_hclust
#'
#' Function to combine smaller clusters from a fast hierarchical algorithm to maximise the BAPS likelihood.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param h a hclust object representing the hierarchical clustering that is to be cut
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' d <- snp_dist(sparse.data)
#' d <- as.dist(d/max(d))
#' h <- hclust(d, method="ward.D2")
#' partition <- best_baps_partition_hclust(sparse.data, h)
#'
#' @export
best_baps_partition_hclust <- function(sparse.data, h){

  # Check inputs
  if(!is.list(sparse.data)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!(class(sparse.data$snp.matrix)=="dgCMatrix")) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.numeric(sparse.data$consensus)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(!is.matrix(sparse.data$prior)) stop("Invalid value for sparse.data! Did you use the import_fasta_sparse_nt function?")
  if(class(h)!="hclust") stop("Invalid value for phylo! Should be a hclust object!")

  llks <- fastbaps:::tree_llk(sparse.data, h$merge)
  n.isolates <- ncol(sparse.data$snp.matrix)

  temp.phylo <- fastbaps:::as.phylo.hclust.node.attributes(h, llks$rk[(n.isolates+1):length(llks$rk)])
  temp.phylo$node.labels <- temp.phylo$node.attributes
  plot(temp.phylo, show.node.label = TRUE)

  node.rk <- c(rep(Inf, n.isolates), temp.phylo$node.attributes)

  cut.nodes <- unlist(lapply(1:(n.isolates+temp.phylo$Nnode), function(n) {
    (node.rk[phangorn::Ancestors(temp.phylo, n, type = "parent")] < 0) && (node.rk[n]>0)
  }))
  cut.nodes <- c(1:(n.isolates+temp.phylo$Nnode))[cut.nodes]

  clusters <- rep(NA, n.isolates)
  clades <- lapply(cut.nodes, function(n) unlist(phangorn:::Descendants(temp.phylo, n, type = "tips")))
  for(i in seq_along(clades)){
    clusters[clades[[i]]] <- i
  }

  return(clusters)
}

# adpated from ape as.phylo
as.phylo.hclust.node.attributes <- function(x, attribute)
{
  N <- dim(x$merge)[1]
  edge <- matrix(0L, 2*N, 2)
  edge.length <- numeric(2*N)
  ## `node' gives the number of the node for the i-th row of x$merge
  node <- integer(N)
  node[N] <- N + 2L
  node.attributes <- rep(NA, N)
  cur.nod <- N + 3L
  j <- 1L
  for (i in N:1) {
    edge[j:(j + 1), 1] <- node[i]
    for (l in 1:2) {
      k <- j + l - 1L
      y <- x$merge[i, l]
      if (y > 0) {
        edge[k, 2] <- node[y] <- cur.nod
        cur.nod <- cur.nod + 1L
        edge.length[k] <- x$height[i] - x$height[y]
        node.attributes[edge[k, 1]-(N+1)] <- attribute[i]
      } else {
        edge[k, 2] <- -y
        edge.length[k] <- x$height[i]
        node.attributes[edge[k, 1]-(N+1)] <- attribute[i]
      }
    }
    j <- j + 2L
  }
  if (is.null(x$labels))
    x$labels <- as.character(1:(N + 1))
  obj <- list(edge = edge, edge.length = edge.length / 2,
              tip.label = x$labels, Nnode = N, node.attributes=node.attributes)
  class(obj) <- "phylo"
  reorder(obj)
}
