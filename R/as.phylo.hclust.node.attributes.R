#' as.phylo.hclust.node.attributes
#'
#' Function to combine smaller clusters from a fast hierarchical algorithm to maximise the BAPS likelihood.
#'
#' @import ape
#'
#' @param x
#' @param attribute
#'
#'
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
