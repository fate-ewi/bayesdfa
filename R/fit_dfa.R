#' Fit a Bayesian DFA
#'
#' @param y A matrix of data to fit. See `data_shape` option to specify whether
#'   this is long or wide format data.
#' @param covar A matrix of covariates, defaults to NULL (not included)
#' @param covar_index A matrix, dimensioned as the number of time series x
#'   number of covariates that indexes which elements of the covariate matrix
#'   are shared across time series. Defaults to a matrix with unique
#'   coefficients estimated for each covariate-time series combination. Elements
#'   may be shared across time series or covariates.
#' @param num_trends Number of trends to fit.
#' @param varIndx Indices indicating which timeseries should have shared
#'   variances.
#' @param zscore Logical. Should the data be standardized first? If not it is
#'   just centered. Centering is necessary because no intercept is included.
#' @param iter Number of iterations in Stan sampling, defaults to 2000.
#' @param thin Thinning rate in Stan sampling, defaults to 1.
#' @param chains Number of chains in Stan sampling, defaults to 4.
#' @param control A list of options to pass to Stan sampling. Defaults to
#'   `list(adapt_delta = 0.99, max_treedepth = 20)`.
#' @param nu_fixed Student t degrees of freedom parameter. If specified as
#'   greater than 100, a normal random walk is used instead of a random walk
#'   with a t-distribution. Defaults to `101`.
#' @param est_correlation Boolean, whether to estimate correlation of
#'   observation error matrix `R`. Defaults to `FALSE`.
#' @param estimate_nu Logical. Estimate the student t degrees of freedom
#'   parameter? Defaults to `FALSE`,
#' @param estimate_trend_ar Logical. Estimate AR(1) parameters on DFA trends?
#'   Defaults to `FALSE``, in which case AR(1) parameters are set to 1
#' @param estimate_trend_ma Logical. Estimate MA(1) parameters on DFA trends?
#'   Defaults to `FALSE``, in which case MA(1) parameters are set to 0.
#' @param sample Logical. Should the model be sampled from? If `FALSE`, then the
#'   data list object that would have been passed to Stan is returned instead.
#'   This is useful for debugging and simulation. Defaults to `TRUE`.
#' @param data_shape If `wide` (the current default) then the input data should
#'   have rows representing the various timeseries and columns representing the
#'   values through time. This matches the MARSS input data format. If `long`
#'   then the input data should have columns representing the various timeseries
#'   and rows representing the values through time.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#' @details Note that there is nothing restricting the loadings and trends from
#'   being inverted (i.e. multiplied by `-1`) for a given chain. Therefore, if
#'   you fit multiple chains, the package will attempt to determine which chains
#'   need to be inverted using the function [find_inverted_chains()].
#' @seealso plot_loadings plot_trends rotate_trends find_swans
#'
#' @export
#'
#' @importFrom rstan sampling
#' @import Rcpp
#' @importFrom graphics lines par plot points polygon segments
#' @importFrom stats na.omit
#'
#' @examples
#' set.seed(42)
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' # only 1 chain and 1000 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, iter = 1000, chains = 1)

fit_dfa <- function(y = y,
                    covar = NULL,
                    covar_index = NULL,
                    num_trends = 1,
                    varIndx = NULL,
                    zscore = TRUE,
                    iter = 2000,
                    chains = 4,
                    thin = 1,
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    nu_fixed = 101,
                    est_correlation = FALSE,
                    estimate_nu = FALSE,
                    estimate_trend_ar = FALSE,
                    estimate_trend_ma = FALSE,
                    sample = TRUE,
                    data_shape = c("wide", "long"),
                    seed = sample.int(.Machine$integer.max, 1),...) {
  data_shape <- match.arg(data_shape)
  if (ncol(y) > nrow(y) && data_shape == "long") {
    warning(
      "ncol(y) > nrow(y) and data_shape == 'long'; are you sure your",
      "input data is in long format?"
    )
  }
  if (ncol(y) < nrow(y) && data_shape == "wide") {
    warning(
      "ncol(y) < nrow(y) and data_shape == 'wide'; are you sure your",
      "input data is in wide format?"
    )
  }
  if (data_shape == "long") {
    y <- t(y)
  }

  if (nrow(y) < 3) {
    stop("fit_dfa() only works with 3 or more time series. We detected ",
      nrow(y), " time series.")
  }

  # parameters for DFA
  N <- ncol(y) # number of time steps
  P <- nrow(y) # number of time series
  K <- num_trends # number of dfa trends
  nZ <- P * K - sum(seq_len(K)) # number of non-zero parameters that are unconstrained

  for (i in seq_len(P)) {
    if (zscore) {
      if (length(unique(na.omit(c(y[i, ])))) == 1L) {
        stop("Can't scale one or more of the time series because all values ",
          "are the same. Remove this/these time series or set `zscore = FALSE`.",
          call. = FALSE
        )
      }
      y[i, ] <- scale(y[i, ], center = TRUE, scale = TRUE)
    } else {
      y[i, ] <- scale(y[i, ], center = TRUE, scale = FALSE)
    }
  }
  Y <- y # included in returned object at end
  # Deal with covariates
  d_covar <- covar

  num_covar <- nrow(d_covar)
  covar_indexing <- covar_index
  if (!is.null(d_covar) && is.null(covar_indexing)) {
    # covariates included but index matrix not, assume independent for all elements
    covar_indexing <- matrix(seq(1, num_covar * P), P, num_covar)
    num_unique_covar <- max(covar_indexing)
  }
  if (is.null(d_covar)) {
    covar_indexing <- matrix(0, P, 0)
    d_covar <- matrix(0, 0, N)
    num_covar <- 0
    num_unique_covar <- 0
  }

  # mat_indx now references the unconstrained values of the Z matrix.
  mat_indx <- matrix(0, P, K)
  start <- 1
  for (k in seq_len(K)) {
    for (p in seq(k + 1, P)) {
      mat_indx[p, k] <- start
      start <- start + 1
    }
  }
  # row_indx and col_indx now references the unconstrained values of the Z matrix.
  row_indx <- matrix((rep(seq_len(P), K)), P, K)[mat_indx > 0]
  col_indx <- matrix(sort(rep(seq_len(K), P)), P, K)[mat_indx > 0]

  diag(mat_indx) <- 1
  row_indx_z <- matrix((rep(seq_len(P), K)), P, K)[mat_indx == 0]
  col_indx_z <- matrix(sort(rep(seq_len(K), P)), P, K)[mat_indx == 0]
  row_indx_z <- c(row_indx_z, 0, 0) # +2 zeros for making stan ok with data types
  col_indx_z <- c(col_indx_z, 0, 0) # +2 zeros for making stan ok with data types
  nZero <- length(row_indx_z)

  # set the model up to have shared variances
  if (is.null(varIndx)) {
    varIndx <- rep(1, P)
  }
  nVariances <- length(unique(varIndx))

  # indices of positive values - Stan can't handle NAs
  row_indx_pos <- matrix(rep(seq_len(P), N), P, N)[!is.na(y)]
  col_indx_pos <- matrix(sort(rep(seq_len(N), P)), P, N)[!is.na(y)]
  n_pos <- length(row_indx_pos)

  row_indx_na <- matrix(rep(seq_len(P), N), P, N)[is.na(y)]
  col_indx_na <- matrix(sort(rep(seq_len(N), P)), P, N)[is.na(y)]
  n_na <- length(row_indx_na)

  y <- y[!is.na(y)]

  # flag for whether to use a normal dist
  use_normal <- if (nu_fixed > 100) 1 else 0
  if (estimate_nu) use_normal <- 0 # competing flags

  # flag for whether to use constraint on Z or not
  # default: only for trends > 1
  zlow = ifelse(num_trends == 1, -100, 0)

  data_list <- list(
    N = N,
    P = P,
    K = K,
    nZ = nZ,
    y = y,
    row_indx = row_indx,
    col_indx = col_indx,
    nZero = nZero,
    varIndx = varIndx,
    nVariances = nVariances,
    row_indx_z = row_indx_z,
    col_indx_z = col_indx_z,
    nZero = nZero,
    row_indx_z = row_indx_z,
    col_indx_z = col_indx_z,
    row_indx_pos = row_indx_pos,
    col_indx_pos = col_indx_pos,
    n_pos = n_pos,
    row_indx_na = row_indx_na,
    col_indx_na = col_indx_na,
    n_na = n_na,
    nu_fixed = nu_fixed,
    d_covar = d_covar,
    num_covar = num_covar,
    covar_indexing = covar_indexing,
    num_unique_covar = num_unique_covar,
    estimate_nu = as.integer(estimate_nu),
    use_normal = use_normal,
    est_cor = as.numeric(est_correlation),
    est_phi = as.numeric(estimate_trend_ar),
    est_theta = as.numeric(estimate_trend_ma),
    zlow = zlow
  )

  pars <- c("x", "Z", "pred", "sigma", "log_lik")
  if (est_correlation) pars <- c(pars, "Omega") # add correlation matrix
  if (!is.null(covar)) pars <- c(pars, "D")
  if (estimate_nu) pars <- c(pars, "nu")
  if (estimate_trend_ar) pars <- c(pars, "phi")
  if (estimate_trend_ma) pars <- c(pars, "theta")

  sampling_args <- list(
    object = stanmodels$dfa,
    data = data_list,
    pars = pars,
    control = control,
    chains = chains,
    iter = iter,
    thin = thin,
    seed = seed,
    ...
  )

  if (sample) {
    mod <- do.call(sampling, sampling_args)
    if (chains > 1) {
      out <- invert_chains(mod, trends = num_trends, print = FALSE)
    } else {
      e <- rstan::extract(mod, permuted = FALSE)
      ep <- rstan::extract(mod, permuted = TRUE)
      out <- list(
        model = mod, samples_permuted = ep, samples = e,
        monitor = rstan::monitor(e)
      )
    }

    out[["data"]] <- Y # keep data included
    out <- structure(out, class = "bayesdfa")
  } else {
    out <- data_list
  }
  out
}

