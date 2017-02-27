#' Summarize Rhat convergence statistics across parameters
#'
#' Pass in rstanfit model object, and optional threshold value for convergence. Returns boolean
#'
#' @param fitted_model A matrix of data to fit. Columns represent time element.
#' @param threshold Threshold for maximum Rhat, default = 1.05
#' @export
#'
#' @importFrom rstan
converge_rhat = function(fitted_model, threshold = 1.05) {
  Rhats = summary(fitted_model)$summary[,"Rhat"]
  return(ifelse(max(Rhats, na.rm=T) < threshold, TRUE, FALSE))
}