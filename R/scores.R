#' Scoring function generator (beta-uniform model)
#' @export
#' @param pvals a numeric vector. P-values that will be used for learning
#'              scoring function
#' @param fdr false discovery rate. It is used for calculating threshold p-value.
#'            P-values less than threshold value will be scored positively.
#' @param threshold_pval custom p-value threshold
#' @param type
#' @param bum_plot whether or not plot BUM fitting plot
#' @param verbose be verbose
#' @param n potential number of genes
#' @seealso \code{\link{posterior_probabilities}}
bum_score <- function(pvals, fdr = 0.01, threshold_pval,
                      type = c("aggressive", "light"), bum_plot = FALSE,
                      n = length(pvals), verbose = FALSE) {
  if(!(class(pvals) == "numeric" && all(pvals >= 0) && all(pvals <= 1))) {
    stop("Invalid p-values")
  }
  type <- match.arg(type, c("aggressive", "light"))

  fit_pvals <- pvals
  if (n > length(pvals)) {
    fit_pvals <- c(fit_pvals, runif(n - length(pvals)))
  }
  fb <- BioNet::fitBumModel(fit_pvals, plot = bum_plot)

  prob <- function(x, fb) {
    sn <- ((1 - fb$lambda) * stats::dbeta(x, fb$a, 1)) / fb$lambda
    sn / (sn + 1)
  }

  if (verbose) {
    print(paste("lambda: ", fb$lambda))
    print(paste("alpha: ", fb$a))
  }
  if (!missing(threshold_pval)) {
    threshold <- threshold_pval
  } else {
    threshold <- BioNet::fdrThreshold(fdr, fb)
  }
  res <- list(p = prob, fb = fb)
  res$score <- if (type == "aggressive") {
    function(x) {
      (fb$a - 1) * (log(x) - log(threshold))
    }
  } else {
    function(x) {
      log(prob(x, fb)) - log(prob(threshold, fb))
    }
  }
  res
}
