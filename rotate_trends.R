rotate_trends = function(fitted_model) {
  
  # Illustrate how to get the trends out of the model
  # get the inverse of the rotation matrix
  Zest = apply(extract(fitted_model)$Z, c(2,3), mean)
  H.inv = varimax(Zest)$rotmat
  
  # rotate factor loadings
  Z.rot = Zest %*% H.inv
  # rotate trends
  states = apply(extract(fitted_model)$x, c(2,3), mean)
  trends.rot = solve(H.inv) %*% states
  
  return(list("Z_rot"=Z.rot, "trends"=t(trends.rot)))
}