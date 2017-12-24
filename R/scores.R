#' @export
BUM_score <- function(pvals, FDR = 0.05, threshold_pval,
                      type = c("aggressive", "light"), bum_plot = FALSE,
                      verbose = FALSE) {
  type <- match.arg(type, c("aggressive", "light"))
  fb <- BioNet::fitBumModel(pvals, plot = bum_plot)
  if (verbose) {
    print(fb)
  }
  if (!missing(threshold_pval)) {
    threshold <- threshold_pval
  } else {
    threshold <- BioNet::fdrThreshold(FDR, fb)
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
