#' Fit multiple models with differing numbers of regimes to trend data
#'
#' @param y Data, time series or trend from fitted DFA model.
#' @param sds Optional time series of standard deviations of estimates. If
#'   passed in, residual variance not estimated.
#' @param min_regimes Smallest of regimes to evaluate, defaults to 1.
#' @param max_regimes Biggest of regimes to evaluate, defaults to 3.
#' @param ... Other parameters to pass to [rstan::sampling()].
#' @param iter MCMC iterations, defaults to 2000.
#' @param thin MCMC thinning rate, defaults to 1.
#' @param chains MCMC chains; defaults to 1 (note that running multiple chains
#'   may result in a "label switching" problem where the regimes are identified
#'   with different IDs across chains).
#' @export
#'
#' @examples
#' data(Nile)
#' find_regimes(log(Nile), iter = 50, chains = 1, max_regimes = 2)

find_regimes <- function(y,
  sds = NULL,
  min_regimes = 1,
  max_regimes = 3,
  iter = 2000,
  thin = 1,
  chains = 1,
  ...) {

  df <- data.frame(regimes = seq(min_regimes, max_regimes), looic = NA)
  best_loo <- 1.0e10
  best_model <- NA
  for (regime in seq(min_regimes, max_regimes)) {
    fit <- fit_regimes(
      y = y, sds = sds, n_regimes = regime, iter = iter, thin = thin,
      chains = chains, ...
    )
    looic <- loo.bayesdfa(fit)
    loo_bad <- loo::pareto_k_table(looic)["(0.7, 1]","Count"]
    loo_very_bad <- loo::pareto_k_table(looic)["(1, Inf)","Count"]
    df$looic[which(df$regimes == regime)] = looic$estimates["looic", "Estimate"]

    if (fit$looic < best_loo) {
      best_loo <- fit$looic
      best_model <- fit
      n_loo_bad <- loo_bad
      n_loo_very_bad <- loo_very_bad
    }
  }

  list(table = df, best_model = best_model, n_loo_bad = n_loo_bad,
    n_loo_very_bad = n_loo_very_bad)
}
