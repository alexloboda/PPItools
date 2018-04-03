extract_component <- function(net, sg) {
  genes <- V(sg)$name
  membership <- components(net)$membership
  id <- which(V(net)$name == genes[1])
  cid <- membership[id]
  remain <- which(membership == cid)
  induced_subgraph(net, remain)
}

#' plots heat-mapped module
#' @export
#' @param g igraph network
#' @param pps named vector of some posterior probabilies
#' @param log_space considering posteriors in log space
#' @param colors number of colors in heat map palette
plot_mwcs <- function(g, pps, log_space = TRUE, colors = 16) {
  if (log_space) {
    pps <- log(pps)
    pps[pps == -Inf] <- min(pps[pps > -Inf])
  }
  cls <- heat.colors(colors)
  pps_cs <- cls[cut(c(pps[V(g)$name], 0), colors)]
  pps_cs <- head(pps_cs, length(pps_cs) - 1)
  V(g)$color <- pps_cs
  V(g)$size <- 8
  V(g)$frame.color <- "white"
  plot(g, layout=layout_with_fr, vertex.label.cex = 0.6)
}

perm_fun <- function(permutator, tf) {
  permutator <- match.arg(permutator, c('igraph', 'ppitools'))
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

#' finds active module of given network
#' @export
#' @inheritParams permutation_freq
#' @param ... other parameters of mwcsr::rmwcs
estimate_network <- function(network,
                             solver_time_factor = 1.0,
                             ...) {
  stopifnot("score" %in% igraph::list.vertex.attributes(network))
  stopifnot(solver_time_factor > 0)

  args <- list(...)
  args$max_iterations <- round(1000 * solver_time_factor)

  solve <- do.call(mwcsr::rmwcs, args)
  solve(network)
}

#' computes frequencies of genes to be in a random network
#' @param pvals named vector of pvals or data frame containing columns \
#' gene, p_value and locus
#' @param network igraph undirected graph with supplied gene names as vertex names
#' @param scf scoring function(e.g. bum_score output)
#' @param permute_interactome should interactome be shuffled
#' @param permute_pvals should p-values be shuffled
#' @param n_permutations number of permutations to do
#' @param permutation_method igraph or ppitools, choose second choice if \
#' you want to change permutation_time_factor
#' @param permutation_time_factor permutator will take time proportionally to this parameter
#' @param solver_time_factor solver will take time proportionally to this parameter
#' @param plot whether to plot a module
#' @param simplify if FALSE sets p-value 0.5 to genes with unknown p-value
#' @param threads number of threads that could be used for computations
#' @export
permutation_freq <- function(pvals,
                             network,
                             scoring_function = bum_score(pvals),
                             permute_interactome = TRUE,
                             permute_pvals = FALSE,
                             n_permutations = 1000,
                             permutation_method = c('igraph', 'ppitools'),
                             permutation_time_factor = 1.0,
                             solver_time_factor = 1.0,
                             plot = FALSE,
                             simplify = TRUE,
                             threads = parallel::detectCores()) {
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  permute <- perm_fun(permutation_method, permutation_time_factor)

  data <- preprocess(pvals, network, simplify, scoring_function)
  lapply(names(data), function(n) assign(n, data[[n]], envir = parent.frame(n = 2)))

  sg <- estimate_network(net, solver_time_factor)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_component(net, sg)


  fs <- foreach::`%dopar%`(
    foreach::foreach(i = 1:n_permutations, .combine = `+`, .packages = "igraph",
                     .inorder = FALSE, .export = "estimate_network"), {
      tmp <- net

      if (permute_pvals) {
        V(tmp)$score <- sample(V(tmp)$score)
      }
      if (permute_interactome) {
        tmp <- permute(tmp)
      }
      module <- estimate_network(tmp, solver_time_factor)
      V(net)$name %in% V(module)$name
    }
  )
  setNames(fs / n_permutations, V(net)$name)
}

#' preprocesses network and p-values
#'
#' It is convinient to do before passing
#' network to estimate_network function, especially in case when set of genes
#' presetnted in pvals does not completely match to set of genes in network
#' @export
#' @inheritParams permutation_freq
preprocess <- function(pvals, network, simplify, scf) {
  if (class(pvals) == "data.frame") {
    pvals <- extract_gwas(pvals)
  }

  if(!(class(pvals) == "numeric" && all(pvals >= 0) && all(pvals <= 1))) {
    stop("Invalid p-values")
  }

  if (length(names(pvals)) == 0) {
    stop("Please name the p-values with gene names ")
  }

  stopifnot(length(pvals) == length(names(pvals)))
  if (!simplify) {
    remain <- setdiff(V(network)$name, names(pvals))
    pvals[remain] <- 0.5;
  }
  network <- simplify(network)

  genes <- intersect(names(pvals), V(network)$name)
  network <- induced_subgraph(network, genes)
  pvals <- pvals[V(network)$name]

  V(network)$score <- unname(scf$score(pvals))

  list(net = network, pvals = pvals, p = scf$p, fb = scf$fb, score = scf$score)
}

#' computes p(U|BM)
#' @import igraph
#' @export
#' @inheritParams permutation_freq
#' @param verbose be verbose
posterior_probabilities <- function(pvals,
                                    network,
                                    scoring_function = BUM_score(pvals),
                                    threads = parallel::detectCores(),
                                    solver_time_factor = 1.0,
                                    simplify = TRUE,
                                    verbose = FALSE,
                                    plot = TRUE) {
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  solve <- mwcsr::rmwcs(max_iterations = round(1000 * solver_time_factor),
                        verbose = verbose)
  data <- preprocess(pvals, network, simplify, scoring_function)
  lapply(names(data), function(n) assign(n, data[[n]], envir = parent.frame(n = 2)))

  sg <- solve(net)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_component(net, sg)

  ns <- V(sg)$name
  sg_weight <- sum(V(sg)$score)

  pps <- foreach::`%dopar%` (
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
  if (plot) plot_mwcs(sg, pps)
  pps
}

#' @export
#' @rdname permutation_freq
#' @param verbose verbosity
BM_pvalue <- function(pvals,
                      network,
                      scoring_function = BUM_score(pvals),
                      threads = parallel::detectCores(),
                      solver_time_factor = 1.0,
                      simplify = TRUE,
                      n_permutations = 100,
                      permutation_method = c('igraph', 'ppitools'),
                      permutation_time_factor = 1.0,
                      verbose = FALSE) {
  permute <- perm_fun(match.arg(permutation_method, c("igraph", "ppitools")),
                      permutation_time_factor)
  cl <- parallel::makeCluster(threads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  solve <- mwcsr::rmwcs(max_iterations = round(1000 * solver_time_factor),
                        verbose = verbose)

  data <- preprocess(pvals, network, simplify, scoring_function)
  lapply(names(data), function(n) assign(n, data[[n]], envir = parent.frame(n = 2)))

  sg <- solve(net)
  if (length(V(sg)) == 0) {
    return(NULL)
  }
  net <- extract_component(net, sg)

  init <- data.frame(weight = c(), degree = c())
  obs <- foreach::`%dopar%` (
    foreach::foreach(i = 1:n_permutations, .combine = rbind, .init = init,
                     .packages = "igraph"), {
      random_net <- permute(net)
      mod <- solve(random_net)
      list(weight = sum(V(mod)$score), degree = length(E(mod)) / length(V(mod)))
  })
  list(p_val_weight = sum(obs$weight >= sum(V(sg)$score)) / permutations,
       p_val_degree = sum(obs$degree >= length(E(sg)) / length(V(sg))) / permutations)
}
