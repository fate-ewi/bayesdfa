#' Fit a Bayesian DFA
#'
#' @param y A matrix of data to fit. Columns represent time element.
#' @param covar A matrix of covariates
#' @param covar_index Indices assigning ??
#' @param num_trends Number of trends to fit.
#' @param varIndx Indices indicating which timeseries should have shared variances.
#' @param zscore Logical. Should the data be standardized first?
#' @param iter Number of iterations in Stan sampling.
#' @param chains Number of chains in Stan sampling.
#' @param control A list of options to pass to Stan sampling.
#' @param nu_fixed Student t degrees of freedom parameter
#' @param tau A fixed parameter describing the standard deviation on the random
#'   walk for the factor loadings in the case of time varying DFA.
#' @param timevarying Logical. If \code{TRUE}, a time varying DFA. Note that the
#'   time varying DFA has not been extensively tested and may not return
#'   sensible answers.
#' @param estimate_nu Logical. Estimate the student t degrees of freedom parameter?
#'
#' @export
#'
#' @importFrom rstan sampling
#' @import Rcpp

fit_dfa = function(y = y,
  covar=NULL,
  covar_index=NULL,
  num_trends = 2,
  varIndx = NULL,
  zscore = TRUE,
  iter = 4000,
  chains = 1,
  control = list(adapt_delta = 0.99),
  nu_fixed = 7,
  tau = 0.1,
  timevarying = FALSE,
  estimate_nu = FALSE) {
  # parameters for DFA
  N = ncol(y)
  P = nrow(y)
  K = num_trends # number of dfa trends
  nZ = P * K - sum(1:K) + K

  if (zscore) {
    for (i in seq_len(P)) {
      y[i, ] <- scale(y[i, ], center = TRUE, scale = TRUE)
    }
  }

  # Deal with covariates
  d_covar = covar

  num_covar = nrow(d_covar)
  covar_indexing = covar_index
  if (!is.null(d_covar) & is.null(covar_indexing)) {
    # covariates included but index matrix not, assume independent for all elements
    covar_indexing = matrix(seq(1, num_covar * P), P, num_covar)
    num_unique_covar = max(covar_indexing)
  }
  if (is.null(d_covar)) {
    covar_indexing = matrix(0, P, 0)
    d_covar = matrix(0, 0, N)
    num_covar = 0
    num_unique_covar = 0
  }

  mat_indx = matrix(0, P, K)
  start = 1
  for (k in 1:K) {
    if (k == 1)
      mat_indx[, k] = (seq_len(nZ))[start:(start + P - k)]
    if (k > 1)
      mat_indx[-c(0:(k - 1)), k] = (seq_len(nZ))[start:(start + P - k)]
    start = start + (P - k + 1)
  }

  row_indx = matrix((rep(seq_len(P), K)), P, K)[which(mat_indx > 0)]
  col_indx = rep(1:K, times = P:(P - K + 1))
  row_indx_z = matrix((rep(seq_len(P), K)), P, K)[which(mat_indx == 0)]
  col_indx_z = matrix(sort(rep(seq_len(K), P)), P, K)[which(mat_indx == 0)]
  row_indx_z = c(row_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  col_indx_z = c(col_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  nZero = length(row_indx_z)

  # set the model up to have shared variances between first two time series,
  # third is different
  if (is.null(varIndx))
    varIndx = rep(1, P)
  nVariances = length(unique(varIndx))

  # indices of positive values - stan can't handle NAs
  row_indx_pos = matrix((rep(seq_len(P), N)), P, N)[which(!is.na(y))]
  col_indx_pos = matrix(sort(rep(seq_len(N), P)), P, N)[which(!is.na(y))]
  n_pos = length(row_indx_pos)
  y = y[which(!is.na(y))]

  data_list = list(
    N,
    P,
    K,
    nZ,
    y,
    row_indx,
    col_indx,
    nZero,
    row_indx_z,
    col_indx_z,
    nZero,
    row_indx_z,
    col_indx_z,
    row_indx_pos,
    col_indx_pos,
    n_pos,
    nu_fixed,
    tau,
    d_covar,
    num_covar,
    covar_indexing,
    num_unique_covar,
    est_df = as.integer(estimate_nu)
  )
  pars <- c("x", "Z", "sigma", "log_lik")
  if (!is.null(covar)) pars <- c(pars, "D")
  if (estimate_nu) pars <- c(pars, "nu")

  if (timevarying) {
    m <- stanmodels$tvdfa_fixed
  } else {
    m <- stanmodels$dfa
  }

  sampling_args <- list(
    object = m,
    data = data_list,
    pars = pars,
    control = control,
    chains = chains,
    iter = iter)

  mod <- do.call(sampling, sampling_args)

  mod
}
