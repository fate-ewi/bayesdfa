#' Estimate the correlation between a DFA trend and some other timeseries
#'
#' Fully incorporates the uncertainty from the posterior of the DFA trend
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}
#' @param y A numeric vector to correlate with the DFA trend. Must be the same
#'   length as the DFA trend.
#' @param trend A number corresponding to which trend to use, defaults to 1.
#' @param time_window Indices indicating a time window slice to use in the
#'   correlation. Defaults to using the entire time window. Can be used to walk
#'   through the timeseries and test the cross correlations.
#' @param trend_samples The number of samples from the trend posterior to use. A
#'   model will be run for each trend sample so this value shouldn't be too
#'   large. Defaults to 100.
#' @param stan_iter The number of samples from the posterior with each Stan
#'   model run, defaults to 300.
#' @param stan_chains The number of chains for each Stan model run, defaults to 1.
#' @param ... Other arguments to pass to \code{\link[rstan]{sampling}}
#'
#' @return A numeric vector of samples from the correlation coefficient posterior.
#'
#' @examples
#' \dontrun{
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' m <- fit_dfa(y = y, num_trends = 1, iter = 600, nu_fixed = 100)
#' r <- rotate_trends(m)
#' correlation <- trend_cor(r, y = rnorm(ncol(r$trends[,1,]), 0, 1))
#' hist(correlation)
#' correlation <- trend_cor(r, y = rnorm(ncol(r$trends[,1,]), 0, 1),
#'   time_window = 5:15)
#' hist(correlation)
#' }
#' @export

trend_cor <- function(rotated_modelfit, y, trend = 1,
  time_window = seq_len(length(y)), trend_samples = 100,
  stan_iter = 300, stan_chains = 1, ...) {

  # must be even to cleanly divide by 2 later:
  if (!stan_iter %% 2 == 0) stan_iter <- stan_iter + 1
  if (max(time_window) > length(y)) stop("Maximum time window value is too large")

  y <- as.numeric(scale(y[time_window]))
  x <- rotated_modelfit$trends[ , trend, time_window] # samples x trend x time

  if (!identical(ncol(x), length(y))) stop("DFA trend and y must be same length")

  samples <- sample(seq_len(nrow(x)), size = trend_samples)
  out <- vapply(seq_len(length(samples)), FUN = function(i) {
    xi <- as.numeric(scale(as.numeric(x[samples[i],])))
    m <- rstan::sampling(stanmodels$corr,
      data = list(x = xi, y = y, N = length(y)),
      iter = stan_iter, chains = stan_chains, warmup = stan_iter/2, ...)
    rstan::extract(m, pars = "beta")[["beta"]]
  }, FUN.VALUE = numeric(length = stan_iter/2))
  out <- as.numeric(out)
  out
}
