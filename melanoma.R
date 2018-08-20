library(igraph)
library(mwcsr)
devtools::load_all()

pvals <- read.table("pvals.txt", header = TRUE, stringsAsFactors = FALSE)

solver <- rmwcs(verbose = TRUE)

interactome_edgelist <- as.matrix(read.table("data/inwebIM_ppi.txt"))
interactome <- graph_from_edgelist(interactome_edgelist, directed = FALSE)
interactome <- interactome - vertex("UBC")
interactome <- simplify(interactome)
interactome <- induced_subgraph(interactome, which(components(interactome)$membership == 1))

colnames(pvals) <- c("gene", "p_value", "locus")
net<-preprocess(pvals=pvals,
                network=interactome,
                simplify = T,
                scf=bum_score(pvals,
                              threshold_pval=5e-8,
                              type="aggressive",
                              bum_plot=FALSE))

instance <- mwcs_instance(net$net, parse_vertex_weights = FALSE)
vertex_weights(instance) <- net$score(net$pvals)

instance <- solve_mwcsp(solver, instance)
plot(solution(instance))

loci <- list()

pvals <- pvals[!is.na(pvals$locus), ]
for (locus in unique(pvals$locus)) {
  if (!is.na(locus)) {
    loci <- c(loci, list(pvals[pvals$locus == locus, ]$gene))
  }
}

tol <- 1e-6
score <- sum(net$score(net$pvals[V(instance$solution)$name]))

all_best <- list()
i <- 1
for (locus in loci) {
   min_score <- score
   all_scores <- c()
   best <- c()
   for (gene in locus) {
     if (gene %in% V(instance$solution)$name) {
       curr_int <- net$net - vertex(gene)
       curr_instance <- mwcs_instance(curr_int, parse_vertex_weights = FALSE)
       curr_pvals <- net$pvals[V(curr_int)$name]
       vertex_weights(curr_instance) <- net$score(curr_pvals)
       curr_instance <- solve_mwcsp(solver, curr_instance)
       curr_score <- sum(net$score(net$pvals[V(curr_instance$solution)$name]))
       all_scores <- c(all_scores, curr_score - score)
       if (curr_score < min_score + tol) {
         best <- c(best, gene)
       }
       if (curr_score < min_score - tol) {
         min_score <- curr_score
         best <- gene
       }
     } else {
       all_scores <- c(all_scores, NaN)
     }
   }
   all_best[[i]] <- list(delta = score - min_score, genes = best, locus = locus,
                         unselected = tail(all_scores[order(all_scores)],
                                                     length(all_scores) - length(best)))
   i <- i + 1
}

for (l in all_best) {
  cat(paste0("Candidate genes: ", do.call(paste, c(list(sep = "/"), as.list(l$locus))),
               "; Selection: ", do.call(paste, c(list(sep = ", "), as.list(l$genes))),
               "; score change: ", -l$delta,
               "; unselected: ", do.call(paste, c(list(sep = ", "), as.list(l$unselected))), "\n"))
}

