extract_component <- function(net, sg) {
  genes <- V(sg)$name
  membership <- components(net)$membership
  id <- which(V(net)$name == genes[1])
  cid <- membership[id]
  remain <- which(membership == cid)
  induced_subgraph(net, remain)
}

binary_search <- function(f, t, tol) {
  l <- 0.0
  r <- 1.0
  while (r - l > tol) {
    m <- (r + l) / 2
    if (f(m) < t) {
      r <- m
    } else {
      l <- m
    }
  }
  (r + l) / 2
}

#' computes posterior smth
#' @import igraph
#' @export
#' @param pvals a named vector of p-values
#' @param network an igraph undirected graph object with vertex attribute "name"
#' @param threads calculations will be run in this number of threads
#' @param verbose be verbose
#' @param simplify if FALSE missing p-values will be filled by 0.5
#' @param solver_time_factor solver of MWCS problem will work amount of time
#'                           proportionally to this factor
#' @param permutation_time_factor time factor for permutations
posterior_probabilities <- function(pvals,
                                    network,
                                    scoring_function = BUM_score(pvals),
                                    threads = parallel::detectCores(),
                                    verbose = FALSE,
                                    simplify = TRUE,
                                    solver_time_factor = 1.0) {
  if(!(class(pvals) == "numeric" && all(pvals >= 0) && all(pvals <= 1))) {
    stop("Invalid p-values")
  }

  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  solve <- mwcsr::rmwcs(max_iterations = round(1000 * solver_time_factor),
                        verbose = verbose)
  network <- simplify(network)

  if (length(names(pvals)) == 0) {
    stop("Please name the p-values with gene names ")
  }
  stopifnot(length(pvals) == length(names(pvals)))

  if (!simplify) {
    remain <- setdiff(V(net)$name, names(pvals))
    pvals[remain] <- 0.5;
  }
  genes <- intersect(names(pvals), V(network)$name)
  net <- induced_subgraph(network, genes)
  pvals <- pvals[V(net)$name]
  score <- scoring_function$score
  p <- scoring_function$p
  fb <- scoring_function$fb
  V(net)$score <- unname(score(pvals))

  sg <- solve(net)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_component(net, sg)

  ns <- V(sg)$name
  sg_weight <- sum(V(sg)$score)

 foreach::`%dopar%` (
    foreach::foreach(gene = ns, .combine = c, .inorder = FALSE,
                     .packages = "igraph", .export = "binary_search"), {
      g <- delete.vertices(net, gene)
      sol <- solve(g)
      diff <- sg_weight - sum(V(sol)$score)
      threshold_weight <- vertex.attributes(net, gene)$score - diff
      threshold_pval <- binary_search(score, threshold_weight, 1e-12)
      p_beta <- stats::pbeta(threshold_pval, fb$a, 1)
      l <- fb$lambda
      p_bum <- l * threshold_pval + (1 - l) * p_beta
      rm(g, sol)
      setNames((1 - p(pvals[gene], fb)) * (threshold_pval / p_bum), gene)
  })
}
