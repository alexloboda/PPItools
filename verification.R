library(igraph)
library(mcmcRanking)
library(stats)
library(mwcsr)
library(parallel)
library(doParallel)

args <- commandArgs(trailingOnly = TRUE)
args <- 3

devtools::load_all(".")

time <- function(a) {
  before <- Sys.time()
  a
  Sys.time() - before
}

sample_module <- function(interactome, n) {
  true_module_vs <- mcmcRanking::sample_subgraph(interactome, module_size = n, niter = 500000)
  induced_subgraph(interactome, true_module_vs)
}

chr_to_int <- function(chr) {
  tmp <- substr(chr, 4, 5)
  if (tmp == "X") {
    return(24)
  }
  if (tmp == "Y") {
    return(25)
  }
  as.integer(tmp)
}

create_locus <- function(gene_locations, gene_name, size) {
  chr <- gene_locations[gene_locations$name2 == gene_name, ]$chrom
  gl <- gene_locations[gene_locations$chrom == chr, ]
  gl <- gl[order(gl$cdsStart), ]
  pos <- which(gl$name2 == gene_name)
  if (length(pos) == 0) {
    return(NULL)
  }
  opt <- pos
  for (i in 1:(size - 1)) {
        opt <- c(pos - i, opt)
  }
  opts <- list()
  for (i in 0:(size - 1)) {
    curr <- opt + i
    if (all(curr >= 1) && all(curr <= nrow(gl))) {
      opts <- c(opts, list(curr))
    }
  }
  s <- sample(opts, 1)
  gl[unlist(s), ]$name2
}

gene_locations <- read.table("data/uniq_genes.txt", header = TRUE, stringsAsFactors = FALSE)
interactome_edgelist <- as.matrix(read.table("data/inwebIM_ppi.txt"))
interactome <- graph_from_edgelist(interactome_edgelist, directed = FALSE)
interactome <- interactome - vertex("UBC")
interactome <- simplify(interactome)
interactome <- induced_subgraph(interactome, which(components(interactome)$membership == 1))

gene_locations <- gene_locations[gene_locations$name2 %in% V(interactome)$name, ]
gene_loc <- data.frame(SNP = gene_locations$name2, BP = (gene_locations$cdsEnd + gene_locations$cdsStart) / 2,
                       CHR = sapply(gene_locations$chrom, chr_to_int))

N <- length(V(interactome))
M <- length(E(interactome))

print(paste(N, M))

size <- as.integer(args[1])
solver <- rmwcs()
alpha <- 0.2
tol = 1e-6

cl <- makeCluster(4)
registerDoParallel(cl)

#count_right <- foreach(i = 1:20, .combine = c, .inorder = FALSE,
#                       .packages = c("igraph", "mwcsr"),
#                       .export = c("bum_score")) %dopar% {
for (i in 1:1) {
    failed <- FALSE
    module <- sample_module(interactome, 200)
    module_genes <- V(module)$name
    instance <- mwcs_instance(interactome, parse_vertex_weights = FALSE)
    pvals <- setNames(runif(N), V(interactome)$name)
    pvals[module_genes] <- rbeta(length(module_genes), alpha, 1.0)
    try <- 0
    while (TRUE) {
      BUM <- bum_score(pvals)
      if (BUM$fb$a < 0.3) {
        break
      } else {
        try <- try + 1
      }
      if (try > 5) {
        return(c())
      }
    }
    loci = list()
    tracking_genes <- names(head(pvals[module_genes][order(pvals[module_genes])], 10))
    print(pvals[tracking_genes])
    for (gene in tracking_genes) {
      locus <- create_locus(gene_locations, gene, size)
      if (is.null(locus)) {
        failed <- TRUE
        break
      }
      pvals[locus] <- pvals[gene]
      loci <- c(loci, list(locus))
    }

    curr_right <- c()

    if (!failed) {
      vertex_weights(instance) <- BUM$score(pvals)
      instance <- solve_mwcsp(solver, instance)
      score <- sum(BUM$score(pvals[V(instance$solution)$name]))
      curr_right <- 0
      for (locus in loci) {
        min_score <- score
        best <- c()
        for (gene in locus) {
          if (gene %in% V(instance$solution)$name) {
            curr_int <- interactome - vertex(gene)
            curr_instance <- mwcs_instance(curr_int, parse_vertex_weights = FALSE)
            curr_pvals <- pvals[V(curr_int)$name]
            vertex_weights(curr_instance) <- BUM$score(curr_pvals)
            curr_instance <- solve_mwcsp(solver, curr_instance)
            curr_score <- sum(BUM$score(pvals[V(curr_instance$solution)$name]))
            if (curr_score < min_score + tol) {
              best <- c(best, gene)
            }
            if (curr_score < min_score - tol) {
              min_score <- curr_score
              best <- gene
            }
          }
        }
        if (length(best) != 0) {
          curr_right <- curr_right + (sum(best %in% tracking_genes) / length(best))
        } else {
          curr_right <- curr_right + 1 / length(locus)
        }
      }
    }
    curr_right
}

stopCluster(cl)

write(count_right, paste0("res", size))
print(count_right)
