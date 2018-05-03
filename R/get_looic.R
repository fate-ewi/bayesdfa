#' Get LOOIC
#'
#' Extract the LOOIC (leave-one-out information criterion) using
#' [loo::loo()].
#'
#' @param fitted_model Output from [fit_dfa()].
#'
#' @export
get_looic <- function(fitted_model) {
  log_lik <- loo::extract_log_lik(fitted_model$model, merge_chains = FALSE)
  rel_eff <- loo::relative_eff(exp(log_lik))
  loo::loo(log_lik, r_eff = rel_eff)$estimates["looic", 1]
}
