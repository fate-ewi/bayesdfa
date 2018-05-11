#' Return LOO information criteria
#'
#' Extract the LOOIC (leave-one-out information criterion) using
#' [loo::loo()]. Note that we've implemented slightly different variants
#' of loo, based on whether the DFA observation model includes correlation
#' between time series or not (default is no correlation). Importantly,
#' these different versions are not directly comparable to evaluate data support
#' for including correlation or not in a DFA. If time series are not correlated,
#' the point-wise log-likelihood for each observation is calculated and used
#' in the loo calculations. However if time series are correlated, then each
#' row of the dataframe (time slice) is assumed to be a joint observation of
#' all variables, and the point-wise log-likelihood is calculated as the
#' joint likelihood of all variables under the multivariate normal distribution.
#'
#' @param fitted_model Output from [fit_dfa()].
#' @param cores Number of cores to use for parallelization.
#'
#' @export
#' @examples
#' \donttest{
#' set.seed(1)
#' s <- sim_dfa(num_trends = 1, num_years = 50, num_ts = 3)
#' m <- fit_dfa(y = s$y_sim, iter = 300, chains = 1, num_trends = 1)
#' loo(m)
#' }

loo.bayesdfa <- function(fitted_model, cores = getOption("mc.cores", 1)) {
  log_lik <- loo::extract_log_lik(fitted_model$model, merge_chains = FALSE)
  rel_eff <- loo::relative_eff(exp(log_lik), cores = cores)
  loo::loo.array(log_lik,
    r_eff = rel_eff,
    cores = cores,
    save_psis = FALSE)
}

#' @name loo
#' @rdname loo.bayesdfa
#' @export
#' @importFrom loo loo
NULL