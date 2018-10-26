#' find_max_cluster_subset
#'
#' Function to find the best partition given a weighting on the nodes
#'
#' @param phylo a phylo object
#' @param node.weights node weights
#'
#' @import Matrix
#'
find_max_cluster_subset <- function(phylo, node.weights){
  n.isolates <- ape::Ntip(phylo)

  nnodes <- phylo$Nnode+n.isolates
  node.index <- 1:nnodes
  node.heights <- ape::node.depth(phylo)

  #initialise
  deltas <- rep(1, nnodes)
  SC.hat <- rep(NA, nnodes)

  #set SC for all leaf nodes
  SC.hat[node.heights==1] <- node.weights[node.heights==1]

  internal.nodes <- node.index[order(node.heights, decreasing = FALSE)]
  internal.node.heights <- node.heights[order(node.heights, decreasing = FALSE)]
  internal.nodes <- internal.nodes[internal.node.heights>1]

  for (n in internal.nodes){
    sum.children <- sum(SC.hat[phangorn::Children(phylo, n)])
    stopifnot(!is.na(sum.children))
    if (node.weights[[n]] <= sum.children){
      SC.hat[[n]] <- sum.children
      deltas[[n]] <- 0
    } else {
      SC.hat[[n]] <- node.weights[[n]]
      deltas[phangorn::Descendants(phylo, n, type = "all")] <- 0
    }
  }

  clusters <- rep(NA, ape::Ntip(phylo))
  names(clusters) <- phylo$tip.label

  i <- 1
  for (n in node.index[deltas==1]){
    clusters[unlist(phangorn::Descendants(phylo, n, type="tips"))] <- i
    i <- i+1
  }

  return(clusters)

}
