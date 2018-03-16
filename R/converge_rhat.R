#' Summarize Rhat convergence statistics across parameters
#'
#' Pass in `rstanfit` model object, and a threshold Rhat value for
#' convergence. Returns boolean.
#'
#' @param fitted_model Samples extracted (with `permuted = FALSE`) from a Stan
#'   model. E.g. output from [invert_chains()].
#' @param threshold Threshold for maximum Rhat.
#' @export
#'
is_converged <- function(fitted_model, threshold = 1.05) {
  Rhats <- fitted_model$monitor[, "Rhat"]
  max(Rhats, na.rm = TRUE) < threshold
}
