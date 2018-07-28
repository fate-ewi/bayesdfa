#' Summarize Rhat convergence statistics across parameters
#'
#' Pass in `rstanfit` model object, and a threshold Rhat value for
#' convergence. Returns boolean.
#'
#' @param fitted_model Samples extracted (with `permuted = FALSE`) from a Stan
#'   model. E.g. output from [invert_chains()].
#' @param threshold Threshold for maximum Rhat.
#' @param parameters Vector of parameters to be included in convergence determination. Defaults = c("sigma","x","Z"). Other elements can be added including "pred", "log_lik", or "lp__"
#' @export
#'
is_converged <- function(fitted_model,
  threshold = 1.05,
  parameters = c("sigma", "x", "Z")) {

  Rhats <-
    fitted_model$monitor[which(grepl(
      paste(parameters, collapse = "|"),
      rownames(fitted_model$monitor)
    ) == TRUE), "Rhat"]

  max(Rhats, na.rm = TRUE) < threshold
}
