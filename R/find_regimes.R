# some code copied from https://github.com/luisdamiano/stancon18
# under CC-BY 4.0

hmm_init <- function(K, x_t) {
  clasif <- kmeans(x_t, K)
  init.mu <- by(x_t, clasif$cluster, mean)
  init.sigma <- by(x_t, clasif$cluster, sd)
  init.order <- order(init.mu)
  list(mu_k = init.mu[init.order], sigma_k = init.sigma[init.order])
}

#' Fit models with differing numbers of regimes to trend data
#'
#' @param y Data, time series or trend from fitted DFA model.
#' @param sds Optional time series of standard deviations of estimates. If passed in, residual variance not estimated
#' @param n_regimes Number of regimes to evaluate
#' @param ... Other parameters to pass to \code{\link[rstan]{sampling}}
#' @param iter MCMC iterations, defaults to 2000
#' @param chains MCMC chains, defaults to 1 (note that running multiple chains
#'   may result in a label switching problem where the regimes are identified
#'   with different IDs across chains).
#' @export
#'
#' @importFrom rstan sampling
#' @import Rcpp
#'
#' @examples
#' \dontrun{
#' data(Nile)
#' find_regimes(log(Nile))
#' }

find_regimes <- function(y, sds = NULL, n_regimes = 2, iter = 2000, chains = 1, ...) {

  est_sigma = 0
  if(is.null(sds)) {
    # estimate sigma, instead of using fixed values
    sds = rep(0, length(y))
    est_sigma = 1
  }

  stan_data = list(
    T = length(y),
    K = n_regimes,
    x_t = y,
    sigma_t = sds,
    est_sigma = est_sigma,
    pars = c("p_1k", "A_ij", "mu_k", "sigma_k", "log_lik", "unalpha_tk", "gamma_tk",
      "unbeta_tk","ungamma_tk","alpha_tk", "beta_tk", "zstar_t", "logp_zstar_t"))

  m <- rstan::sampling(stanmodels$hmm_gaussian,
    data = stan_data,
    iter = iter,
    chains = chains,
    init = function() {hmm_init(n_regimes, y)},
    ...)

  list(model = m, y = y)
}

#' Plot the state probabilities from \code{find_regimes}
#'
#' @param model A model returned by \code{\link{find_regimes}}.
#' @param probs A numeric vector of quantiles to plot the credible intervals at.
#' @param regime_prob_threshold The probability density that must be above 0.5
#'   before we classify a regime (only affects \code{"means"} plot).
#' @details
#' Note that the original timeseries data (dots) are shown scaled between 0 and 1.
#'
#' @export
#' @examples
#' \dontrun{
#' data(Nile)
#' m <- find_regimes(log(Nile))
#' plot_regime_model(m)
#' plot_regime_model(m, type = "means")
#'
#' set.seed(1)
#' y <- c(rnorm(20, 0, 0.2), rnorm(20, 0.5, 0.2), rnorm(20, 1, 0.2))
#' m <- find_regimes(y, n_regimes = 2)
#' plot_regime_model(m)
#' plot_regime_model(m, type = "means", regime_prob_threshold = 0.95)
#'
#' set.seed(1)
#' y <- c(rnorm(20, 0, 0.1), rnorm(20, 0.5, 0.1), rnorm(20, 1, 0.1))
#' m <- find_regimes(y, n_regimes = 3)
#' plot_regime_model(m)
#' plot_regime_model(m, type = "means", regime_prob_threshold = 0.95)
#' }
plot_regime_model <- function(model, probs = c(0.05, 0.95),
  type = c("probability", "means"),
  regime_prob_threshold = 0.9) {

  gamma_tk <- rstan::extract(model$model, pars = 'gamma_tk')[[1]]
  mu_k <- rstan::extract(model$model, pars = 'mu_k')[[1]]
  l <- apply(gamma_tk, 2:3, quantile, probs = probs[[1]])
  u <- apply(gamma_tk, 2:3, quantile, probs = probs[[2]])
  med <- apply(gamma_tk, 2:3, quantile, probs = 0.5)
  range01 <- function(x) (x-min(x))/(max(x)-min(x))
  mu_k_low <- apply(mu_k, 2, quantile, probs = probs[[1]])
  mu_k_high <- apply(mu_k, 2, quantile, probs = probs[[2]])
  mu_k <- apply(mu_k, 2, median)
  confident_regimes <- apply(gamma_tk, 2:3, function(x) mean(x > 0.5) > regime_prob_threshold)
  regime_indexes <- apply(confident_regimes, 1, function(x) {
    w <- which(x)
    ifelse(length(w) == 0, NA, w)
  })

  if (type[[1]] == "probability") {
    oldpar <- par("mfrow")
    par(mfrow = c(1, ncol(med)))
    for (i in seq_len(ncol(med))) {
      plot(l[,i], ylim = c(0, 1), col = "grey40", lty = 2, type = "n",
        main = paste("State", LETTERS[i]), ylab = "Probability of being in given state",
        xlab = "Time")
      polygon(c(1:nrow(u), nrow(u):1), c(l[,i], rev(u[,i])), col = "grey70", border = "grey70")
      lines(1:nrow(u), med[,i], col = "black", lwd = 2)
      points(1:nrow(u), range01(model$y), col = "#FF000070", pch = 3)

    }
    par(mfrow = oldpar)
  } else {
    plot(as.numeric(model$y), col = "#FF000070", pch = 3, ylab = "Time series value",
      xlab = "Time")
    if (!all(is.na(regime_indexes))) {
      for (i in seq_along(regime_indexes)) {
        segments(x0 = i-0.5, x1 = i+0.5, y0 = mu_k[regime_indexes[i]], y1 = mu_k[regime_indexes[i]])
        polygon(c(i-0.5, i-0.5, i+0.5, i+0.5),
          c(mu_k_low[regime_indexes[i]], mu_k_high[regime_indexes[i]],
            mu_k_high[regime_indexes[i]], mu_k_low[regime_indexes[i]]),
          border = NA, col = "#00000050")
      }
    }
  }
}
