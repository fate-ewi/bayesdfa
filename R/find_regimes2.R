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
#' @param n_regimes Number of regimes to evaluate
#' @param ... Other parameters to pass to \code{\link[rstan]{sampling}}
#' @param iter MCMC iterations, defaults to 2000
#' @param chains MCMC chains, defaults to 1 (note that running multiple chains
#'   may result in a label switching problem where the regimes are identified
#'   with different IDs across chains).
#' @export
#'
#' @examples
#' \dontrun{
#' data(Nile)
#' find_regimes2(log(Nile))
#' }

find_regimes2 <- function(y, n_regimes = 2, iter = 2000, chains = 1, ...) {

  stan_data = list(
    T = length(y),
    K = n_regimes,
    x_t = y)

  m <- rstan::sampling(stanmodels$hmm_gaussian,
    data = stan_data,
    iter = iter,
    chains = chains,
    init = function() {hmm_init(n_regimes, y)},
    ...)

  list(model = m, y = y)
}

#' Plot the state probabilities from \code{find_regimes2}
#'
#' @param model A model returned by \code{\link{find_regimes2}}.
#' @param probs A numeric vector of quantiles to plot the credible intervals at.
#'
#' @details
#' Note that the original timeseries data (dots) are shown scaled between 0 and 1.
#'
#' @export
#' @examples
#' \dontrun{
#' data(Nile)
#' m <- find_regimes2(log(Nile))
#' plot_regime_model(m)
#' }
plot_regime_model <- function(model, probs = c(0.1, 0.9)) {

  gamma_tk <- rstan::extract(model$model, pars = 'gamma_tk')[[1]]
  l <- apply(gamma_tk, 2:3, quantile, probs = probs[[1]])
  u <- apply(gamma_tk, 2:3, quantile, probs = probs[[2]])
  med <- apply(gamma_tk, 2:3, quantile, probs = 0.5)
  range01 <- function(x) (x-min(x))/(max(x)-min(x))

  par(mfrow = c(1, ncol(med)))
  for (i in seq_len(ncol(med))) {
    plot(l[,i], ylim = c(0, 1), col = "grey40", lty = 2, type = "n",
      main = paste("State", LETTERS[i]), ylab = "Probability of being in given state",
      xlab = "Time")
    polygon(c(1:nrow(u), nrow(u):1), c(l[,i], rev(u[,i])), col = "grey70", border = "grey70")
    lines(1:nrow(u), med[,i], col = "black", lwd = 2)
    points(1:nrow(u), range01(model$y), col = "#FF000070", pch = 3)
  }
  par(mfrow = c(1, 1))
}
