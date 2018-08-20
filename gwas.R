library(igraph)
library(mwcsr)
devtools::load_all()

rm(list = ls())

time <- function(x) {
  before <- Sys.time()
  force(x)
  after <- Sys.time()
  after - before
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

gene_locations <- read.table("data/uniq_genes.txt", header = TRUE, stringsAsFactors = FALSE)
gene_locations <- gene_locations[gene_locations$name2 %in% V(interactome)$name, ]
gene_loc <- data.frame(SNP = gene_locations$name2, BP = (gene_locations$cdsEnd + gene_locations$cdsStart) / 2,
                       CHR = sapply(gene_locations$chrom, chr_to_int))

build_distribution <- function(priors, gene_candidates, true_genes, add = FALSE,
                               color = "darkgrey", title = "") {
  for (id in unique(gene_candidates$snp_id)) {
      ns <- gene_candidates[gene_candidates$snp_id == id,]$gene
      ns <- intersect(ns, names(priors))
      if (sum(priors[ns]) == 0.0) {
        priors[ns] <- 1 / length(ns)
      } else {
        priors[ns] <- priors[ns] / sum(priors[ns])
      }
  }
  ps <- priors[true_genes]
  ps <- ps[!is.na(ps)]
  mean <- sum(ps)
  d <- sum(sapply(ps, function(x) x * (1 - x)))
  max <- length(ps)
  t <- setNames(sapply(0:max, function(x) stats::dnorm(x, mean, d)), 0:max)
  barplot(t, main = title, add = add, col = color, ylim = c(0, 0.5), xlim = c(0, 20), xlab = "true positives",
          ylab = "probability")
}

filter_genes <- function(df, set, evidence) {
  snps <- df[df$gene %in% set, ]$snp_id
  replaced <- df[df$gene %in% set, ]
  replaced$evidence <- evidence
  rbind(replaced, df[!(df$snp_id %in% snps), ])
}

# fitting
t <- read.table("UC.PC5.assoc.logistic", header = TRUE)
t <- t[!is.na(t$P),]
#t <- t[t$pval < 0.01 | ((2 * t$AC) / t$nCompleteSamples) > 0.05, ]
fit_pvals <- t$P
fb <- bum_score(fit_pvals, bum_plot = TRUE)

# building network
net_df <- read.table("data/inwebIM_ppi.txt")[, c(1, 2)]
net <- igraph::graph_from_edgelist(as.matrix(net_df))

# reading labeled nearest & GRB genes
df <- read.table("nat_all_and_gws_annot.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")
df <- df[df$gene != "", ]
df[df$pval < 1e-50, ]$pval <- 1e-50

df$p_value <- df$pval

# generates data frame containing genes, not snps
# GRB has higher priority than nearest
# coding genes not considered(only as part of a locus)
# head(genes) for more info
# do not remove i <- 1, it is sufficient part(bad design tho, who cares)
i <- 1
genes <- plyr::adply(df, .margins = 1, .fun = function(x) {
  genes_str <- if (is.na(x$grb_candidate_genes)) x$nearest_gene_candidate_genes
               else x$grb_candidate_genes
  gs <- strsplit(genes_str, ";")[[1]]
  res <- data.frame(gene = gs, stringsAsFactors = FALSE)
  res$evidence <- if(is.na(x$grb_candidate_genes)) "nearest" else "GRB"
  res$p_value <- x$p_value
  i <<- i + 1
  res
}, .expand = FALSE, .id = "locus")

# removing duplicated entries
evid_imp <- setNames(2:1, c("GRB", "nearest_gene"))
genes <- genes[order(evid_imp[genes$evidence], genes$p, decreasing = TRUE), ]
genes <- genes[!duplicated(genes$gene),]

# reading mark and nature true positives
mark <- scan("gold", what = character())
nature <- scan("TP", what = character())
union_true <- union(mark, nature)
locuses <- genes[genes$gene %in% union_true,]$locus
dataset <- genes[genes$locus %in% locuses, ]


##############

solver <- rmwcs()

interactome_edgelist <- as.matrix(read.table("data/inwebIM_ppi.txt"))
interactome <- graph_from_edgelist(interactome_edgelist, directed = FALSE)
interactome <- interactome - vertex("UBC")
interactome <- simplify(interactome)
interactome <- induced_subgraph(interactome, which(components(interactome)$membership == 1))

net<-preprocess(pvals=genes,
                network=interactome,
                simplify = T,
                scf=bum_score(genes,
                              threshold_pval=5e-8,
                              type="aggressive",
                              bum_plot=FALSE))

instance <- mwcs_instance(net$net, parse_vertex_weights = FALSE)
vertex_weights(instance) <- net$score(net$pvals)

instance <- solve_mwcsp(solver, instance)
plot(solution(instance))

loci <- list()
pvals <- genes

pvals <- pvals[!is.na(pvals$locus), ]
for (locus in unique(pvals$locus)) {
  if (locus %in% locuses) {
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


