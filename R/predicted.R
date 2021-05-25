#' Calculate predicted value from DFA object
#'
#' Pass in `rstanfit` model object. Returns array of predictions, dimensioned
#' number of MCMC draws x number of MCMC chains x time series length x number of time series
#'
#' @param fitted_model Samples extracted (with `permuted = FALSE`) from a Stan
#'   model. E.g. output from [invert_chains()].
#' @export
#' @examples
#' \dontrun{
#' set.seed(42)
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' # only 1 chain and 1000 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, iter = 2000, chains = 3, num_trends = 1)
#' pred <- predicted(m)
#' }
predicted <- function(fitted_model) {
  Z <- rstan::extract(fitted_model$model, "Z", permuted = FALSE)
  x <- rstan::extract(fitted_model$model, "x", permuted = FALSE)
  Zperm <- rstan::extract(fitted_model$model, "Z", permuted = TRUE)
  xperm <- rstan::extract(fitted_model$model, "x", permuted = TRUE)

  n_ts <- dim(Zperm$Z)[2]
  n_y <- dim(xperm$x)[3]
  n_chains <- dim(Z)[2]
  n_trends <- dim(xperm$x)[2]
  n_mcmc <- dim(x)[1]

  pred <- array(0, c(n_mcmc, n_chains, n_y, n_ts))
  for (i in 1:n_mcmc) {
    for (chain in 1:n_chains) {
      # for each MCMC draw / chain
      x_i <- t(matrix(x[i, chain, ], nrow = n_trends, ncol = n_y))
      Z_i <- t(matrix(Z[i, chain, ], nrow = n_ts, ncol = n_trends))
      pred[i, chain, , ] <- x_i %*% Z_i
    }
  }

  return(pred)
}
