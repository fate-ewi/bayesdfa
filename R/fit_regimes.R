# Some of the following code copied (and modified) from
# https://github.com/luisdamiano/stancon18
# under CC-BY 4.0

#' Fit models with differing numbers of regimes to trend data
#'
#' @param y Data, time series or trend from fitted DFA model.
#' @param sds Optional time series of standard deviations of estimates.
#'   If passed in, residual variance not estimated. Defaults to `NULL`.
#' @param n_regimes Number of regimes to evaluate, defaults 2
#' @param ... Other parameters to pass to [rstan::sampling()].
#' @param iter MCMC iterations, defaults to 2000.
#' @param thin MCMC thinning rate, defaults to 1.
#' @param chains MCMC chains, defaults to 1 (note that running multiple chains
#'   may result in a label switching problem where the regimes are identified
#'   with different IDs across chains).
#' @export
#'
#' @importFrom rstan sampling
#' @importFrom loo extract_log_lik loo
#' @import Rcpp
#'
#' @examples
#' data(Nile)
#' fit_regimes(log(Nile), iter = 50, n_regimes = 1)

fit_regimes <- function(y,
  sds = NULL,
  n_regimes = 2,
  iter = 2000,
  thin = 1,
  chains = 1,
  ...) {

  est_sigma <- 0
  if (is.null(sds)) {
    # estimate sigma, instead of using fixed values
    sds <- rep(0, length(y))
    est_sigma <- 1
  }

  if (n_regimes < 1) stop("`n_regimes` must be an integer >= 1.", call. = FALSE)

  if (identical(as.integer(n_regimes), 1L)) {
    stan_data <- list(
      T = length(y),
      K = 1,
      x_t = y,
      sigma_t = sds,
      est_sigma = est_sigma,
      pars = c("mu_k", "sigma_k", "log_lik")
    )

    m <- rstan::sampling(object=stanmodels$regime_1,
      data = stan_data,
      iter = iter,
      chains = chains,
      init = function() {
        hmm_init(n_regimes, y)
      },
      ...
    )
  }

  if (n_regimes > 1) {
    stan_data <- list(
      T = length(y),
      K = n_regimes,
      x_t = y,
      sigma_t = sds,
      est_sigma = est_sigma,
      pars = c(
        "p_1k", "A_ij", "mu_k", "sigma_k", "log_lik", "unalpha_tk", "gamma_tk",
        "unbeta_tk", "ungamma_tk", "alpha_tk", "beta_tk", "zstar_t",
        "logp_zstar_t"
      )
    )

    m <- rstan::sampling(object=stanmodels$hmm_gaussian,
      data = stan_data,
      iter = iter,
      thin = thin,
      chains = chains,
      init = function() {
        hmm_init(n_regimes, y)
      },
      ...
    )
  }

  log_lik <- loo::extract_log_lik(m, merge_chains=FALSE)
  #n_chains = dim(rstan::extract(m, "log_lik", permuted=FALSE))[2]
  rel_eff <- loo::relative_eff(exp(log_lik))
  # calculate looic
  looic <- loo::loo(log_lik, r_eff = rel_eff)$estimates["looic",1]

  list(model = m, y = y, looic = looic)
}
