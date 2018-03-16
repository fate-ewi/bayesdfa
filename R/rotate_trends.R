#' Rotate the trends from a DFA
#'
#' @param fitted_model Output from [fit_dfa()].
#' @param conf_level Probability level for CI.
#' @importFrom stats median quantile sd
#'
#' @export
#'
#' @importFrom rstan extract
#' @examples
#' \dontrun{
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' m <- fit_dfa(y = y, num_trends = 1, iter = 600, chains = 1)
#' r <- rotate_trends(m)
#' }

rotate_trends = function(fitted_model, conf_level = 0.95) {

  # get the inverse of the rotation matrix
  n_mcmc = dim(fitted_model$samples)[2] * dim(fitted_model$samples)[1]

  temp <- reshape_samples(fitted_model$samples)
  Z <- temp$Z
  x <- temp$x
  n_ts = dim(Z)[2]
  n_trends = dim(x)[2]
  n_years = dim(x)[3]

  # do rotation for each MCMC draw (slow)
  mcmc_trends_rot = array(0, dim = c(n_mcmc, n_trends, n_years))
  mcmc_Z_rot = array(0, dim = c(n_mcmc, n_ts, n_trends))
  if(n_trends > 1) {
    for(i in seq_len(n_mcmc)) {
      Zest = Z[i,,]
      H.inv = varimax(Zest)$rotmat
      # rotate factor loadings
      Z.rot = Zest %*% H.inv
      mcmc_Z_rot[i,,] = Z.rot
      # rotate trends
      states = x[i,,]
      trends.rot = solve(H.inv) %*% states
      mcmc_trends_rot[i,,] = trends.rot
    }
   }
   if(n_trends==1) {
     mcmc_trends_rot = x
     mcmc_Z_rot = Z
   }

    list(
      Z_rot = mcmc_Z_rot,
      trends = mcmc_trends_rot,
      Z_rot_mean = apply(mcmc_Z_rot, c(2, 3), mean),
      Z_rot_median = apply(mcmc_Z_rot, c(2, 3), median),
      trends_mean = apply(mcmc_trends_rot, c(2, 3), mean),
      trends_median = apply(mcmc_trends_rot, c(2, 3), median),
      trends_lower = apply(mcmc_trends_rot, c(2, 3),
        quantile, (1 - conf_level) / 2),
      trends_upper = apply(mcmc_trends_rot, c(2, 3),
        quantile, 1 - (1 - conf_level) / 2)
    )

}

# Reshape samples:
reshape_samples <- function(samp) {
  s <- reshape2::melt(samp)
  z <- dplyr::filter(s, grepl("Z\\[", .data$parameters))
  z$trend <- as.numeric(gsub("Z\\[[0-9]+,([0-9]+)\\]", "\\1", z$parameters))
  z$ts <- as.numeric(gsub("Z\\[([0-9]+),([0-9]+)\\]", "\\1", z$parameters))
  Z <- reshape2::acast(z, iterations + chains ~ ts ~ trend, value.var = "value")
  x <- dplyr::filter(s, grepl("x\\[", .data$parameters))
  x$trend <- as.numeric(gsub("x\\[([0-9]+),([0-9]+)\\]", "\\1", x$parameters))
  x$time <- as.numeric(gsub("x\\[([0-9]+),([0-9]+)\\]", "\\2", x$parameters))
  x <- reshape2::acast(x, iterations + chains ~ trend ~ time, value.var = "value")
  list(Z = Z, x = x)
}
