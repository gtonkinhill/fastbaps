#' best_baps_partition
#'
#' Function to combine smaller clusters from a fast hierarchical algorithm to maximise the BAPS likelihood.
#'
#' @import Matrix
#'
#' @param sparse.data a sparse SNP data object returned from import_fasta_sparse_nt
#' @param phylo a phylo object representing the hierarchical clustering that is to be cut
#'
#' @return a final clustering
#'
#' @examples
#'
#' fasta.file.name <- system.file("extdata", "seqs.fa", package = "fastbaps")
#' sparse.data <- import_fasta_sparse_nt(fasta.file.name)
#' sim.matrix <- snp_similarity(sparse.data)
#' x <- as.dist(1-sim.matrix/max(sim.matrix))
#' phylo <- ape::as.phylo(hclust(x, method="average"))
#' partition <- best_baps_partition(sparse.data, phylo)
#'
#' phylo <- phytools::read.newick(system.file("extdata", "seqs.fa.treefile", package = "fastbaps"))
#' partition <- best_baps_partition(sparse.data, phylo)
#'
#' @export
best_baps_partition <- function(sparse.data, phylo){
  n.isolates <- ape::Ntip(phylo)

  #reorder matrx to match leaves of phylo
  sparse.data$snp.matrix <- sparse.data$snp.matrix[,match(phylo$tip.label, colnames(sparse.data$snp.matrix))]

  node.partitions <- lapply(1:(phylo$Nnode+n.isolates),  function(n) {
    return(unlist(phangorn::Descendants(phylo, n, type="tips")))
  })

  node.mls <- llks <- fastbaps:::part_llks(sparse.data, partitions = node.partitions)$llk

  clusters <- fastbaps:::find_max_cluster_subset(phylo, node.mls)

  return(clusters)
}
