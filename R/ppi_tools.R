extract_subgraph <- function(net, sg) {
  genes <- V(sg)$name
  membership <- components(net)$membership
  id <- which(V(net)$name == genes[1])
  cid <- membership[id]
  remain <- which(membership == cid)
  induced_subgraph(net, remain)
}

#' @import igraph
#' @export
posterior_probabilities <- function(gene_pvals,
                                    network,
                                    permutations = 10000,
                                    scoring_function = BUM_score(gene_pvals),
                                    threads = parallel::detectCores()) {
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  solve <- mwcsr::rmwcs()

  if (length(names(gene_pvals)) == 0) {
    stop("Please name the p-values with gene names ")
  }
  stopifnot(length(gene_pvals) == length(names(gene_pvals)))

  genes <- intersect(names(gene_pvals), V(network)$name)
  net <- induced_subgraph(network, genes)
  gene_pvals <- gene_pvals[V(net)$name]
  V(net)$score <- unname(scoring_function(gene_pvals))

  sg <- solve(net)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_subgraph(net, sg)

  res <- foreach::`%dopar%`(
    foreach::foreach(i = 1:permutations, .combine = `+`, .inorder = FALSE), {
      random_net <- igraph::sample_degseq(igraph::degree(net), method = "vl")
      igraph::vertex.attributes(random_net) <- igraph::vertex.attributes(net)
      igraph::V(sg)$name %in% igraph::V(solve(random_net))$name
    }
  )

  setNames(res / permutations, V(sg)$name)
}
