#' Rotate the trends from a DFA
#'
#' @param fitted_model Output from \code{\link{fit_dfa}}.
#' @param conf_level Probability level for CI.
#' @importFrom stats median quantile
#'
#' @export
#'
#' @importFrom rstan extract

rotate_trends = function(fitted_model, conf_level = 0.95) {

  # get the inverse of the rotation matrix
  n_mcmc = dim(extract(fitted_model)$Z)[1]
  Z = extract(fitted_model)$Z
  x = extract(fitted_model)$x
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
      trends_lower = apply(mcmc_trends_rot, c(2, 3), quantile, (1 - conf_level) / 2),
      trends_upper = apply(mcmc_trends_rot, c(2, 3), quantile, 1 - (1 - conf_level) / 2)
    )

}
