#' Summarize Rhat convergence statistics across parameters
#'
#' Pass in rstanfit model object, and optional threshold value for
#' convergence. Returns boolean.
#'
#' @param fitted_model Samples extracted (with permuted = FALSE) from a Stan model.
#'   E.g. output from \code{\link{invert_chains}}.
#' @param threshold Threshold for maximum Rhat, default = 1.05
#' @export
#'
is_converged <- function(fitted_model, threshold = 1.05) {
  # Rhats <- rstan::summary(fitted_model)$summary[, "Rhat"]
  Rhats <- fitted_model$monitor[, "Rhat"]
  max(Rhats, na.rm = TRUE) < threshold
}
