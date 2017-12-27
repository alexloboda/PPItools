#' @import igraph
#' @param graph accompanied by edge attribute `class` (1..n)
shake <- function(graph, number_of_permutations = length(E(graph)) * 10,
                  hard_stop = number_of_permutations * 10) {
  edge_list <- as_edgelist(graph, names = FALSE)
  df = data.frame(from = edge_list[, 1],
                  to = edge_list[, 2])
  num_nodes <- length(V(graph))
  res <- shake_internal(df, num_nodes, number_of_permutations, hard_stop)
  g <- graph_from_edgelist(as.matrix(res), directed = FALSE)
  edge.attributes(g) <- edge.attributes(graph)
  vertex.attributes(g) <- vertex.attributes(graph)
  g
}
