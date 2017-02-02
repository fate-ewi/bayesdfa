rotate_trends = function(fitted_model) {
  
  # Illustrate how to get the trends out of the model
  # get the inverse of the rotation matrix
  n_mcmc = dim(extract(fitted_model)$Z)[1]
  Z = extract(fitted_model)$Z
  x = extract(fitted_model)$x
  n_ts = dim(Z)[2]
  n_trends = dim(x)[2]
  n_years = dim(x)[3]
  mcmc_trends_rot = array(0, dim = c(n_mcmc, n_trends, n_years))
  mcmc_Z_rot = array(0, dim = c(n_mcmc, n_ts, n_trends))
  for(i in 1:n_mcmc) {
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
  return(list("Z_rot"=mcmc_Z_rot, "trends"=mcmc_trends_rot))
}