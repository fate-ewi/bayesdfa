#' Create initial values for the HMM model.
#'
#' @param K The number of regimes or clusters to fit. Called by [rstan::sampling()].
#' @param x_t A matrix of values. Called by [rstan::sampling()].
#' @importFrom stats kmeans
#'
#' @return list of initial values (mu, sigma)
hmm_init <- function(K, x_t) {
  clasif <- kmeans(x_t, K)
  init.mu <- by(x_t, clasif$cluster, mean)
  init.sigma <- by(x_t, clasif$cluster, sd)
  init.order <- order(init.mu)
  list(mu_k = init.mu[init.order], sigma_k = init.sigma[init.order])
}
