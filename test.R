rm(list = ls())
devtools::load_all()
library(igraph)

time <- function(x) {
  before <- Sys.time()
  force(x)
  after <- Sys.time()
  after - before
}

filter_genes <- function(df, set, evidence) {
  snps <- df[df$gene %in% set, ]$snp_id
  replaced <- df[df$gene %in% set, ]
  replaced$evidence <- evidence
  rbind(replaced, df[!(df$snp_id %in% snps), ])
}

net_df <- read.table("scores")[, c(1, 2)]
net <- igraph::graph_from_edgelist(as.matrix(net_df), directed = FALSE)
#
#fit_pvals <- read.table("UC.PC5.assoc.logistic", header = TRUE)
#fit_pvals <- fit_pvals$P
#fit_pvals <- fit_pvals[!is.na(fit_pvals)]
pvals <- read.table("fsgs_pvals")
pvals <- setNames(pvals[[2]], pvals[[1]])

#score <- bum_score(fit_pvals, type = "aggressive", bum_plot = T, threshold_pval = 5e-8)
score <- bum_score(pvals, type = "aggressive", bum_plot = T, fdr = 0.1)

#df <- read.table("nat_all_and_gws_annot.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#df <- df[df$gene != "", ]
#df[df$pval < 1e-50, ]$pval <- 1e-50
#
#i <- 1
#genes <- plyr::adply(df, .margins = 1, .fun = function(x) {
#  genes_str <- if (is.na(x$grb_candidate_genes)) x$nearest_gene_candidate_genes
#               else x$grb_candidate_genes
#  gs <- strsplit(genes_str, ";")[[1]]
#  res <- data.frame(gene = gs, stringsAsFactors = FALSE)
#  res$evidence <- if(is.na(x$grb_candidate_genes)) "nearest" else "GRB"
#  res$p <- x$p
#  i <<- i + 1
#  res
#}, .expand = FALSE, .id = "snp_id")
#
#mark <- scan("gold", what = character())
#nature <- scan("TP", what = character())
#
#genes <- filter_genes(genes, mark, "mark")
#genes <- filter_genes(genes, nature, "nature")
#
#genes <- genes[genes$evidence != "nearest", ]
#pvals <- setNames(genes$p, genes$gene)

print(BM_pvalue(pvals, net, score, simplify = TRUE, threads = 4))
pps <- posterior_probabilities(pvals, net, score, simplify = TRUE, threads = 4)

for (gene in dimnames(pps)[[2]]) {
  png(paste0("plots_melanoma/", gene))
  t <- pps[, gene]
  last <- max(which(t != 0))
  data <- setNames(t[1:last], 0:(last - 1))
  barplot(data, main = gene)
  dev.off()
}

