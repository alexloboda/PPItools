#' Scoring function generator (beta-uniform model)
#' @export
#' @param pvals a numeric vector.Pp-values that will be used for learning
#'              scoring function
#' @param fdr false discovery rate. It is used for calculating threshold p-value.
#'            P-values less than threshold value will be scored positively.
#' @param threshold_pval custom p-value threshold
#' @param type
#' @param bum_plot whether or not plot BUM fitting plot
#' @param verbose be verbose
#' @seealso \code{\link{posterior_probabilities}}
bum_score <- function(pvals, fdr = 0.01, threshold_pval,
                      type = c("aggressive", "light"), bum_plot = FALSE,
                      verbose = FALSE) {
  if(!(class(pvals) == "numeric" && all(pvals >= 0) && all(pvals <= 1))) {
    stop("Invalid p-values")
  }
  type <- match.arg(type, c("aggressive", "light"))
  fb <- BioNet::fitBumModel(pvals, plot = bum_plot)
  if (verbose) {
    print(paste("lambda: ", fb$lambda))
    print(paste("alpha: ", fb$a))
  }
  if (!missing(threshold_pval)) {
    threshold <- threshold_pval
  } else {
    threshold <- BioNet::fdrThreshold(fdr, fb)
  }
  if (type == "aggressive") {
    function(x) {
      (fb$a - 1) * (log(x) - log(threshold))
    }
  } else {
    prob <- function(x, fb) {
      sn <- ((1 - fb$lambda) * stats::dbeta(x, fb$a, 1)) / fb$lambda
      sn / (sn + 1)
    }
    function(x) {
      log(prob(x, fb)) - log(prob(threshold, fb))
    }
  }
}
