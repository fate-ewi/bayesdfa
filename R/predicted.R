#' Calculate predicted value from DFA object
#'
#' Pass in `rstanfit` model object, and a threshold Rhat value for
#' convergence. Returns boolean.
#'
#' @param fitted_model Samples extracted (with `permuted = FALSE`) from a Stan
#'   model. E.g. output from [invert_chains()].
#' @export
#'
predicted <- function(fitted_model) {
  Z = rstan::extract(fitted_model$model, "Z", permuted=FALSE)
  x = rstan::extract(fitted_model$model, "x", permuted=FALSE)

  n_chains = dim(Z)[2]
  n_ts = dim(fitted_model$data)[1]
  n_y = dim(fitted_model$data)[2]
  n_trends = dim(Z)[3]/n_ts
  n_mcmc = dim(x)[1]

  pred = array(0, c(n_mcmc, n_chains, n_y, n_ts))
  for(i in 1:n_mcmc) {
    for(chain in 1:n_chains) {
      # for each MCMC draw / chain
      x_i = t(matrix(x[i,chain,], nrow=n_trends, ncol=n_y))
      Z_i = t(matrix(Z[i,chain,], nrow=n_ts, ncol=n_trends))
      pred[i,chain,,] = x_i %*% Z_i
    }
  }

  return(pred)
}
