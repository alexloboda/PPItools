extract_component <- function(net, sg) {
  genes <- V(sg)$name
  membership <- components(net)$membership
  id <- which(V(net)$name == genes[1])
  cid <- membership[id]
  remain <- which(membership == cid)
  induced_subgraph(net, remain)
}

perm_fun <- function(permutator, tf) {
  if (permutator == "igraph") {
    if (tf != 1.0) {
      stop("igraph permutator doesn't support time control")
    }
    f <- function(g) {
      random_net <- igraph::sample_degseq(igraph::degree(g), method = "vl")
      igraph::vertex.attributes(random_net) <- igraph::vertex.attributes(g)
      random_net
    }
  } else {
    f <- function(g) {
      shake(g, number_of_permutations = round(length(E(g)) * 10 * tf))
    }
  }
  f
}

pvalue_test <- function(pdf, g, ps) {
  m <- nrow(pdf)
  setNames(sapply(V(g)$name, function(name) {
    pdf[(degree(g, name) + 1), name] / ps
  }), V(g)$name)
}

#' computes posterior smth
#' @import igraph
#' @export
#' @param gene_pvals a named vector of p-values
#' @param network an igraph undirected graph object with vertex attribute "name"
#' @param threads calculations will be run in this number of threads
#' @param verbose be verbose
#' @param simplify if FALSE missing p-values will be filled by 0.5
#' @param solver_time_factor solver of MWCS problem will work amount of time
#'                           proportionally to this factor
#' @param permutation_method a method that will be used for permuations. 'package'
#'                           supports different time factors
#' @param permutation_time_factor time factor for permutations
posterior_probabilities <- function(gene_pvals,
                                    network,
                                    permutations = 1000,
                                    scoring_function = BUM_score(gene_pvals),
                                    threads = parallel::detectCores(),
                                    verbose = FALSE,
                                    simplify = TRUE,
                                    solver_time_factor = 1.0,
                                    permutation_method = c("igraph", "package"),
                                    permutations_time_factor = 1.0) {
  pvals <- gene_pvals
  if(!(class(pvals) == "numeric" && all(pvals >= 0) && all(pvals <= 1))) {
    stop("Invalid p-values")
  }

  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  permute <- perm_fun(match.arg(permutation_method, c("igraph", "package")),
                                permutations_time_factor)
  solve <- mwcsr::rmwcs(max_iterations = round(1000 * solver_time_factor),
                        verbose = verbose)
  network <- simplify(network)

  if (length(names(gene_pvals)) == 0) {
    stop("Please name the p-values with gene names ")
  }
  stopifnot(length(gene_pvals) == length(names(gene_pvals)))

  if (!simplify) {
    remain <- setdiff(V(net)$name, names(gene_pvals))
    gene_pvals[remain] <- 0.5;
  }
  genes <- intersect(names(gene_pvals), V(network)$name)
  net <- induced_subgraph(network, genes)
  gene_pvals <- gene_pvals[V(net)$name]
  V(net)$score <- unname(scoring_function(gene_pvals))

  sg <- solve(net)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_component(net, sg)
  degree_pdf <- matrix(0, ncol = length(V(sg)), nrow = length(V(net)),
                       dimnames = list(c(), col = V(sg)$name))

  comb <- function(acc, x) {
    pos <- cbind(x + 1, seq_len(length(x)))
    acc[pos] <- acc[pos] + 1
    acc
  }

  if (permutations == 0) {
    return(setNames(rep(0, length(V(sg))), V(sg)$name))
  }

  res <- foreach::`%dopar%`(
    foreach::foreach(i = 1:permutations, .init = degree_pdf, .combine = comb,
                     .inorder = FALSE), {
      random_net <- permute(net)
      mod <- solve(random_net)
      res <- setNames(rep(0, length(V(sg))), V(sg)$name)
      int <- intersect(V(mod)$name, V(sg)$name)
      res[int] <- degree(mod, int)
      res
    }
  )
  pvalue_test(res, sg, permutations)
}
