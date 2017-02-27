#' Summarize Rhat convergence statistics across parameters
#'
#' Pass in rstanfit model object, and optional threshold value for convergence. Returns boolean
#'
#' @param fitted_model A Stan model. E.g. a model from \code{\link{fit_dfa}}.
#' @param threshold Threshold for maximum Rhat, default = 1.05
#' @export
#'
is_converged <- function(fitted_model, threshold = 1.05) {
  Rhats <- rstan::summary(fitted_model)$summary[ ,"Rhat"]
  max(Rhats, na.rm = TRUE) < threshold
}
