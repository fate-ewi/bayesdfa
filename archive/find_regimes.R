#' Fit models with differing numbers of regimes to trend data
#'
#' Pass in rstanfit model object, and optional threshold value for
#' convergence. Returns boolean.
#'
#' @param y Data, time series or trend from fitted DFA model.
#' @param lower Lower bound of number of regimes to evaluate, defaults to 1
#' @param upper Upper bound of number of regimes to evaluate, defaults to 5
#' @param mcmc_iter MCMC iterations, defaults to 2000
#' @param mcmc_chains MCMC chains, defaults to 3
#' @param threshold Threshold for maximum Rhat, default = 1.05
#' @export
#'
#' @examples
#' \dontrun{
#' data(Nile)
#' find_regimes(y=log(Nile), upper=3)
#' }
#' @importFrom loo loo extract_log_lik
#'
find_regimes = function(y,
  lower = 1,
  upper = 5,
  mcmc_iter = 2000,
  mcmc_chains = 3,
  threshold = 1.05) {
  out_df = data.frame("Model" = 1:upper,
    "LOOIC" = NA,
    "converged" = FALSE)
  n = length(y)
  pars = c("u", "sigma", "log_lik")

  if (lower == 1) {
    # fit null model with 1 regime
    data_list = list(N = n, y = y)
    fitted_model = stan(
      "exec/regime_1.stan",
      data = data_list,
      iter = mcmc_iter,
      chains = mcmc_chains
    )
    out_df$LOOIC[1] = loo::loo(loo::extract_log_lik(fitted_model))$looic
    Rhats = summary(fitted_model)$summary[, "Rhat"]
    out_df$converged[1] = max(Rhats, na.rm = TRUE) < threshold
  }
  if (lower <= 2 & upper >= 2) {
    # fit null model with 2 regimes
    data_list = list(N = n, y = y)
    fitted_model = stan(
      "exec/regime_2.stan",
      data = data_list,
      iter = mcmc_iter,
      chains = mcmc_chains
    )
    out_df$LOOIC[2] = loo::loo(loo::extract_log_lik(fitted_model))$looic
    Rhats = summary(fitted_model)$summary[, "Rhat"]
    out_df$converged[2] = max(Rhats, na.rm = TRUE) < threshold
  }
  if (upper >= 3) {
    lower_bound = ifelse(lower <= 3, 3, lower)
    for (regimes in lower_bound:upper) {
      data_list = list(
        N = n,
        n_regime = regimes,
        y = y,
        ones = rep(1, regimes - 1)
      )
      fitted_model = stan(
        "exec/regime_3plus.stan",
        data = data_list,
        iter = mcmc_iter,
        chains = mcmc_chains
      )
      out_df$LOOIC[regimes] = loo::loo(loo::extract_log_lik(fitted_model))$looic
      Rhats = summary(fitted_model)$summary[, "Rhat"]
      out_df$converged[regimes] = max(Rhats, na.rm = TRUE) < threshold
    }
  }

  return(out_df)
}
