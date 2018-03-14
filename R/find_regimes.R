#' Fit multiple models with differing numbers of regimes to trend data
#'
#' @param y Data, time series or trend from fitted DFA model.
#' @param sds Optional time series of standard deviations of estimates. If
#'   passed in, residual variance not estimated.
#' @param min_regimes Smallest of regimes to evaluate.
#' @param max_regimes Biggest of regimes to evaluate.
#' @param ... Other parameters to pass to [rstan::sampling()].
#' @param iter MCMC iterations.
#' @param chains MCMC chains; defaults to 1 (note that running multiple chains
#'   may result in a "label switching" problem where the regimes are identified
#'   with different IDs across chains).
#' @export
#'
#' @examples
#' \dontrun{
#' data(Nile)
#' find_regimes(log(Nile))
#' }
find_regimes = function(y, sds = NULL, min_regimes = 1, max_regimes = 3,
  iter = 2000, chains = 1, ...) {
  df <- data.frame(regimes = seq(min_regimes, max_regimes), looic = NA)
  best_loo <- 1.0e10
  best_model <- NA
  for (regime in seq(min_regimes, max_regimes)) {
    fit <- fit_regimes(
      y = y, sds = sds, n_regimes = regime, iter = iter,
      chains = chains, ...
    )
    df$looic[which(df$regimes == regime)] <- fit$looic
    if (fit$looic < best_loo) {
      best_loo <- fit$looic
      best_model <- fit
    }
  }

  list(table = df, best_model = best_model)
}
