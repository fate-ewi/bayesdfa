#' Fit a Bayesian DFA
#'
#' @param y A matrix of data to fit. See `data_shape` option to specify whether
#'   this is long or wide format data. Wide format data (default) is a matrix with
#'   time across columns and unique time series across rows, and can only contain 1 observation
#'   per time series - time combination. In contrast, long format data is a data frame that includes
#'   observations ("obs"), time ("time") and time series ("ts") identifiers -- the benefit of long
#'   format is that multiple observations per time series can be included. Correlation matrix currently
#'   not estimated if data shape is long.
#' @param num_trends Number of trends to fit.
#' @param varIndx Indices indicating which timeseries should have shared
#'   variances.
#' @param scale Character string, used to standardized data. Can be "zscore" to center and
#' standardize data, "center" to just standardize data, or "none". Defaults to "zscore"
#' @param iter Number of iterations in Stan sampling, defaults to 2000. Used for both
#' [rstan::sampling()] and [rstan::vb()]
#' @param thin Thinning rate in Stan sampling, defaults to 1.
#' @param chains Number of chains in Stan sampling, defaults to 4.
#' @param control A list of options to pass to Stan sampling. Defaults to
#'   `list(adapt_delta = 0.99, max_treedepth = 20)`.
#' @param nu_fixed Student t degrees of freedom parameter. If specified as
#'   greater than 100, a normal random walk is used instead of a random walk
#'   with a t-distribution. Defaults to `101`.
#' @param est_correlation Boolean, whether to estimate correlation of
#'   observation error matrix `R`. Defaults to `FALSE`. Currently can't be estimated if data are in long format.
#' @param estimate_nu Logical. Estimate the student t degrees of freedom
#'   parameter? Defaults to `FALSE`,
#' @param estimate_trend_ar Logical. Estimate AR(1) parameters on DFA trends?
#'   Defaults to `FALSE``, in which case AR(1) parameters are set to 1
#' @param estimate_trend_ma Logical. Estimate MA(1) parameters on DFA trends?
#'   Defaults to `FALSE``, in which case MA(1) parameters are set to 0.
#' @param estimate_process_sigma Logical. Defaults FALSE, whether or not to estimate process error sigma. If not estimated,
#'   sigma is fixed at 1, like conventional DFAs.
#' @param equal_process_sigma Logical. If process sigma is estimated, whether or not to estimate a single shared value across trends (default)
#'   or estimate equal values for each trend
#' @param estimation Character string. Should the model be sampled using [rstan::sampling()] ("sampling",default),
#' [rstan::optimizing()]("optimizing"), variational inference [rstan::vb()]("vb"),
#' or no estimation done ("none"). No estimation may be useful for debugging and simulation.
#' @param data_shape If `wide` (the current default) then the input data should
#'   have rows representing the various timeseries and columns representing the
#'   values through time. This matches the MARSS input data format. If `long`
#'   then the long format data is a data frame that includes observations ("obs"),
#'   time ("time") and time series ("ts") identifiers -- the benefit of long
#'   format is that multiple observations per time series can be included
#' @param obs_covar Optional dataframe of data with 4 named columns ("time","timeseries","covariate","value"), representing: (1) time, (2) the time series
#'   affected, (3) the covariate number for models with more than one covariate affecting each
#'   trend, and (4) the value of the covariate
#' @param pro_covar Optional dataframe of data with 4 named columns ("time","trend","covariate","value"), representing: (1) time, (2) the trend
#'   affected, (3) the covariate number for models with more than one covariate affecting each
#'   trend, and (4) the value of the covariate
#' @param z_bound Optional hard constraints for estimated factor loadings -- really only applies to model with 1 trend. Passed in as a 2-element vector representing the lower and upper bound, e.g. (0, 100) to constrain positive
#' @param z_model Optional argument allowing for elements of Z to be constrained to be proportions (each time series modeled as a mixture of trends). Arguments can be "dfa" (default) or "proportion"
#' @param trend_model Optional argument to change the model of the underlying latent trend. By default this is set to 'rw', where the trend
#' is modeled as a random walk - as in conentional DFA. Alternative options are 'spline', where B-splines are used to model the trends
#' or 'gp', where gaussian predictive processes are used. If models other than 'rw' are used, there are some key points. First, the MA and AR
#' parameters on these models will be turned off. Second, for B-splines the process_sigma becomes an optional scalar on the spline coefficients,
#' and is turned off by default. Third, the number of knots can be specified (more knots = more wiggliness, and n_knots < N). For models
#' with > 2 trends, each trend has their own spline coefficients estimated though the knot locations are assumed shared. If knots aren't specified,
#' the default is N/3.
#' @param n_knots The number of knots for the B-spline of Gaussian predictive process models. Optional, defaults to round(N/3)
#' @param knot_locs Locations of knots (optional), defaults to uniform spacing between 1 and N
#' @param family String describing the observation model. Default is "gaussian",
#'   but included options are "gamma", "lognormal", negative binomial ("nbinom2"),
#'   "poisson", or "binomial". The binomial family is assumed to have logit link,
#'   gaussian family is assumed to be identity, and the rest are log-link.
#' @param gp_theta_prior A 2-element vector controlling the prior on the Gaussian process parameter in cov_exp_quad.
#'   This prior is a half-Student t prior, with the first argument of gp_theta_prior being the degrees of freedom (nu),
#'   and the second element being the standard deviation
#' @param expansion_prior Defaults to FALSE, if TRUE uses the parameter expansion prior of Ghosh & Dunson 2009
#' @param ... Any other arguments to pass to [rstan::sampling()].
#' @param par_list A vector of parameter names of variables to be estimated by Stan. If NULL, this will default to
#'   c("x", "Z", "sigma", "log_lik", "psi","xstar") for most models -- though if AR / MA, or Student-t models are used
#'   additional parameters will be monitored. If you want to use diagnostic tools in rstan, including moment_matching,
#'   you will need to pass in a larger list. Setting this argument to "all" will monitor all parameters, enabling the use
#'   of diagnostic functions -- but making the models a lot larger for storage. Finally, this argument may be a custom string
#'   of parameters to monitor, e.g. c("x","sigma")
#' @param verbose Whether to print iterations and information from Stan, defaults to FALSE.
#' @details Note that there is nothing restricting the loadings and trends from
#'   being inverted (i.e. multiplied by `-1`) for a given chain. Therefore, if
#'   you fit multiple chains, the package will attempt to determine which chains
#'   need to be inverted using the function [find_inverted_chains()].
#' @seealso plot_loadings plot_trends rotate_trends find_swans
#'
#' @export
#'
#' @importFrom rstan sampling optimizing vb
#' @importFrom splines bs
#' @importFrom stats dist gaussian
#' @import Rcpp
#' @importFrom graphics lines par plot points polygon segments
#' @importFrom stats na.omit runif
#'
#' @examples
#' set.seed(42)
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' # only 1 chain and 250 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, iter = 50, chains = 1)
#' \dontrun{
#' # example of observation error covariates
#' set.seed(42)
#' obs_covar <- expand.grid("time" = 1:20, "timeseries" = 1:3, "covariate" = 1)
#' obs_covar$value <- rnorm(nrow(obs_covar), 0, 0.1)
#' m <- fit_dfa(y = s$y_sim, iter = 50, chains = 1, obs_covar = obs_covar)
#'
#' # example of process error covariates
#' pro_covar <- expand.grid("time" = 1:20, "trend" = 1:2, "covariate" = 1)
#' pro_covar$value <- rnorm(nrow(pro_covar), 0, 0.1)
#' m <- fit_dfa(y = s$y_sim, iter = 50, chains = 1, num_trends = 2, pro_covar = pro_covar)
#'
#' # example of long format data
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' obs <- c(s$y_sim[1, ], s$y_sim[2, ], s$y_sim[3, ])
#' long <- data.frame("obs" = obs, "ts" = sort(rep(1:3, 20)), "time" = rep(1:20, 3))
#' m <- fit_dfa(y = long, data_shape = "long", iter = 50, chains = 1)
#'
#' # example of long format data with obs covariates
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' obs <- c(s$y_sim[1, ], s$y_sim[2, ], s$y_sim[3, ])
#' long <- data.frame("obs" = obs, "ts" = sort(rep(1:3, 20)), "time" = rep(1:20, 3))
#' obs_covar <- expand.grid("time" = 1:20, "timeseries" = 1:3, "covariate" = 1:2)
#' obs_covar$value <- rnorm(nrow(obs_covar), 0, 0.1)
#' m <- fit_dfa(y = long, data_shape = "long", iter = 50, chains = 1, obs_covar = obs_covar)
#'
#' # example of model with Z constrained to be proportions and wide format data
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' m <- fit_dfa(y = s$y_sim, z_model = "proportion", iter = 50, chains = 1)
#'
#' # example of model with Z constrained to be proportions and long format data
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' obs <- c(s$y_sim[1, ], s$y_sim[2, ], s$y_sim[3, ])
#' long <- data.frame("obs" = obs, "ts" = sort(rep(1:3, 20)), "time" = rep(1:20, 3))
#' m <- fit_dfa(y = long, data_shape = "long", z_model = "proportion", iter = 50, chains = 1)
#'
#' #' # example of B-spline model with wide format data
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' m <- fit_dfa(y = s$y_sim, iter = 50, chains = 1, trend_model = "spline", n_knots = 10)
#'
#' # example of B-spline model with wide format data
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' m <- fit_dfa(y = s$y_sim, iter = 50, chains = 1, trend_model = "gp", n_knots = 5)
#' }
fit_dfa <- function(y = y,
                    num_trends = 1,
                    varIndx = NULL,
                    scale = c("zscore", "center", "none"),
                    iter = 2000,
                    chains = 4,
                    thin = 1,
                    control = list(adapt_delta = 0.99, max_treedepth = 20),
                    nu_fixed = 101,
                    est_correlation = FALSE,
                    estimate_nu = FALSE,
                    estimate_trend_ar = FALSE,
                    estimate_trend_ma = FALSE,
                    estimate_process_sigma = FALSE,
                    equal_process_sigma = TRUE,
                    estimation = c("sampling", "optimizing", "vb", "none"),
                    data_shape = c("wide", "long"),
                    obs_covar = NULL,
                    pro_covar = NULL,
                    z_bound = NULL,
                    z_model = c("dfa", "proportion"),
                    trend_model = c("rw", "spline", "gp"),
                    n_knots = NULL,
                    knot_locs = NULL,
                    par_list = NULL,
                    family = "gaussian",
                    verbose = FALSE,
                    gp_theta_prior = c(3, 1),
                    expansion_prior = FALSE,
                    ...) {
  # check arguments
  data_shape <- match.arg(data_shape, c("wide", "long"))
  z_model <- match.arg(z_model, c("dfa", "proportion"))
  trend_model <- match.arg(trend_model, c("rw", "spline", "gp"))

  obs_model <- match(family, c(
    "gaussian", "gamma", "poisson", "nbinom2",
    "binomial", "lognormal"
  ))
  if (is.na(obs_model)) {
    stop("Error: family not found. Please enter family as gaussian(), gamma(), etc.")
  }
  if (family != "gaussian") {
    if (data_shape == "wide") stop("Error: if family is non-gaussian, data must be in long format")
    if (est_correlation == TRUE) stop("Error: correlation can't be estimated with non-gaussian data. Please set est_correlation=FALSE")
  }

  orig_data <- y # save original data

  if (ncol(y) < nrow(y) && data_shape[1] == "wide") {
    warning(
      "ncol(y) < nrow(y) and data_shape == 'wide'; are you sure your",
      "input data is in wide format?"
    )
  }
  if (data_shape[1] == "long") {
    if (est_correlation == TRUE) {
      stop("Estimation of the observation error correlation matrix not currently estimated when data are in long format")
    }
    if (length(which(names(y) == "ts")) == 0) {
      stop("Error: data shape is long, and must contain a field 'ts' representing time series dimension")
    }
    if (length(which(names(y) == "time")) == 0) {
      stop("Error: data shape is long, and must contain a field 'time' representing time dimension")
    }
    if (length(which(names(y) == "obs")) == 0) {
      stop("Error: data shape is long, and must contain a field 'obs' representing observations")
    }
    # rescale if needed
    # y$time <- y$time - min(y[["time"]]) + 1 # min time now = 1
    y$ts <- as.numeric(as.factor(y[["ts"]]))
    N <- max(y[["time"]])
    P <- max(y[["ts"]])
  }

  if (data_shape[1] == "wide") {
    N <- ncol(y) # number of time steps
    P <- nrow(y) # number of time series
    if (nrow(y) < 3) {
      stop(
        "fit_dfa() only works with 3 or more time series. We detected ",
        nrow(y), " time series."
      )
    }
  }

  if (!is.null(obs_covar)) {
    if (ncol(obs_covar) != 4) {
      stop("observation covariates must be in a data frame with 4 columns")
    }
  }
  if (!is.null(pro_covar)) {
    if (ncol(pro_covar) != 4) {
      stop("process covariates must be in a data frame with 4 columns")
    }
  }
  if (!is.null(z_bound) && length(z_bound) != 2) {
    stop("if you're using z bounds, this needs to be a 2-element vector")
  }

  # parameters for DFA
  K <- num_trends # number of dfa trends
  nZ <- P * K - sum(seq_len(K)) # number of non-zero parameters that are unconstrained

  # standardizing data by rows only works if data provided in "wide" format
  if (family == "gaussian") {
    if (data_shape[1] == "wide") {
      for (i in seq_len(P)) {
        if (scale[1] == "zscore") {
          if (length(unique(na.omit(c(y[i, ])))) == 1L) {
            stop("Can't scale one or more of the time series because all values ",
              "are the same. Remove this/these time series or set `scale` = `center`.",
              call. = FALSE
            )
          }
          y[i, ] <- scale(y[i, ], center = TRUE, scale = TRUE)
        }
        if (scale[1] == "center") {
          y[i, ] <- scale(y[i, ], center = TRUE, scale = FALSE)
        }
      }
    } else {
      if (scale[1] == "zscore") {
        # standardize
        for (i in seq_len(P)) {
          indx <- which(y[["ts"]] == i)
          y[indx, "obs"] <- scale(y[indx, "obs"], center = TRUE, scale = TRUE)
        }
      }
      if (scale[1] == "center") {
        # just center
        for (i in seq_len(P)) {
          indx <- which(y[["ts"]] == i)
          y[indx, "obs"] <- scale(y[indx, "obs"], center = TRUE, scale = FALSE)
        }
      }
    }
  }
  Y <- y # included in returned object at end

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

  # row_indx_z and col_indx_z contain locations of zeros in Z matrix of loadings
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
  if (data_shape[1] == "wide") {
    row_indx_pos <- matrix(rep(seq_len(P), N), P, N)[!is.na(y)]
    col_indx_pos <- matrix(sort(rep(seq_len(N), P)), P, N)[!is.na(y)]
    n_pos <- length(row_indx_pos)
    row_indx_na <- matrix(rep(seq_len(P), N), P, N)[is.na(y)]
    col_indx_na <- matrix(sort(rep(seq_len(N), P)), P, N)[is.na(y)]
    n_na <- length(row_indx_na)
    y <- y[!is.na(y)]
  } else {
    y <- y[which(!is.na(y[["obs"]])), ]
    row_indx_pos <- y[["ts"]]
    col_indx_pos <- y[["time"]]
    n_pos <- length(row_indx_pos)
    # these are just dummy placeholders
    row_indx_na <- matrix(1, 1, 1)[is.na(runif(1))]
    col_indx_na <- matrix(1, 1, 1)[is.na(runif(1))]
    n_na <- length(row_indx_na)
    y <- y[["obs"]]
  }

  # flag for whether to use a normal dist
  use_normal <- if (nu_fixed > 100) 1 else 0
  if (estimate_nu) use_normal <- 0 # competing flags

  # covariates
  if (!is.null(obs_covar)) {
    obs_covar_index <- as.matrix(obs_covar[, c("time", "timeseries", "covariate")])
    num_obs_covar <- nrow(obs_covar_index)
    n_obs_covar <- length(unique(obs_covar_index[, "covariate"]))
    obs_covar_value <- obs_covar[, "value"]

    if (data_shape[1] == "wide") {
      match_obs_covar <- rep(0, num_obs_covar)
    } else {
      match_obs_covar <- match(paste(obs_covar$time, obs_covar$timeseries), paste(Y$time[which(!is.na(Y$obs))], Y$ts[which(!is.na(Y$obs))]))
      keep <- which(!is.na(match_obs_covar))
      # keep covariates not associated with missing values
      obs_covar_index <- obs_covar_index[keep, ]
      num_obs_covar <- nrow(obs_covar_index)
      n_obs_covar <- length(unique(obs_covar_index[, "covariate"]))
      obs_covar_value <- obs_covar[keep, "value"]
      # keep matches
      match_obs_covar <- match_obs_covar[keep]
    }
  } else {
    num_obs_covar <- 0
    n_obs_covar <- 0
    obs_covar_value <- c(0)[0]
    match_obs_covar <- c(0)[0]
    obs_covar_index <- matrix(0, 1, 3)[c(0)[0], ]
  }
  if (!is.null(pro_covar)) {
    pro_covar_index <- as.matrix(pro_covar[, c("time", "trend", "covariate")])
    num_pro_covar <- nrow(pro_covar_index)
    n_pro_covar <- length(unique(pro_covar_index[, "covariate"]))
    pro_covar_value <- pro_covar[, "value"]
  } else {
    num_pro_covar <- 0
    n_pro_covar <- 0
    pro_covar_value <- c(0)[0]
    pro_covar_index <- matrix(0, 1, 3)[c(0)[0], ]
  }

  if (is.null(z_bound)) {
    z_bound <- c(-100, 100)
  }

  n_sigma_process <- 1
  if (equal_process_sigma == FALSE) n_sigma_process <- K
  est_sigma_process <- 0
  if (estimate_process_sigma == TRUE) est_sigma_process <- 1

  # default args that need to be passed in
  est_spline <- 0
  est_gp <- 0
  est_rw <- 1 # these are flags specifying model structure. default is rw
  if (is.null(n_knots)) n_knots <- round(N / 3)
  if (is.null(knot_locs)) knot_locs <- seq(1, N, length.out = n_knots)
  distKnots <- matrix(0, n_knots, n_knots)
  # distKnots21 <- matrix(0, N, n_knots)
  distKnots21_pred <- rep(0, n_knots)
  # set up cubic b-splines design matrix
  B_spline <- matrix(0, n_knots, N)

  if (trend_model == "spline") {
    est_spline <- 1
    est_rw <- 0
    # turn of things conventionally estimated when trend is a random walk
    estimate_trend_ar <- FALSE
    estimate_trend_ma <- FALSE
    estimate_nu <- FALSE
    B_spline <- t(splines::bs(1:N, df = n_knots, degree = 3, intercept = TRUE))
  }
  if (trend_model == "gp") {
    # Gaussian kernel
    est_gp <- 1
    est_rw <- 0
    if (is.null(knot_locs)) knot_locs <- seq(1, N, length.out = n_knots)
    distKnots <- as.matrix(stats::dist(knot_locs)) # distances between time stamps
    distAll <- as.matrix(stats::dist(c(1:N, knot_locs))) # distances between data and knot locs
    # distKnots21 <- t(distAll[-seq_len(N), 1:N])
    distKnots21_pred <- as.matrix(stats::dist(c(N + 1, knot_locs)))[1, -1]
    # distKnots <- distKnots ^ 2
    # distKnots21 <- distKnots21 ^ 2
    distKnots21_pred <- distKnots21_pred^2
    est_sigma_process <- 1 # turn this on as a scale for variance
    estimate_trend_ar <- FALSE
    estimate_trend_ma <- FALSE
    estimate_nu <- FALSE
  }

  y_int <- rep(0, length(y))
  if (family %in% c("binomial", "nbinom2", "poisson")) {
    y_int <- as.integer(y)
  }
  est_sigma_params <- ifelse(family %in% c("gaussian", "lognormal"), 1, 0)
  est_gamma_params <- ifelse(family == "gamma", 1, 0)
  est_nb2_params <- ifelse(family == "nbinom2", 1, 0)

  data_list <- list(
    N = N,
    P = P,
    K = K,
    nZ = nZ,
    y = y,
    y_int = y_int,
    row_indx = row_indx,
    col_indx = col_indx,
    nZero = nZero,
    varIndx = varIndx,
    nVariances = nVariances,
    row_indx_z = row_indx_z,
    col_indx_z = col_indx_z,
    row_indx_pos = row_indx_pos,
    col_indx_pos = col_indx_pos,
    n_pos = n_pos,
    row_indx_na = row_indx_na,
    col_indx_na = col_indx_na,
    n_na = n_na,
    nu_fixed = nu_fixed,
    estimate_nu = as.integer(estimate_nu),
    use_normal = use_normal,
    est_cor = as.numeric(est_correlation),
    est_phi = as.numeric(estimate_trend_ar),
    est_theta = as.numeric(estimate_trend_ma),
    num_obs_covar = num_obs_covar,
    n_obs_covar = n_obs_covar,
    obs_covar_value = obs_covar_value,
    obs_covar_index = obs_covar_index,
    match_obs_covar = match_obs_covar,
    num_pro_covar = num_pro_covar,
    n_pro_covar = n_pro_covar,
    pro_covar_value = pro_covar_value,
    pro_covar_index = pro_covar_index,
    z_bound = z_bound,
    long_format = ifelse(data_shape[1] == "wide", 0, 1),
    proportional_model = ifelse(z_model[1] == "dfa", 0, 1),
    est_sigma_process = est_sigma_process,
    n_sigma_process = n_sigma_process,
    est_rw = est_rw,
    est_spline = est_spline,
    B_spline = B_spline,
    n_knots = n_knots,
    knot_locs = knot_locs,
    est_gp = est_gp,
    # distKnots = distKnots,
    # distKnots21 = distKnots21,
    obs_model = obs_model,
    distKnots21_pred = matrix(distKnots21_pred, nrow = 1),
    est_sigma_params = est_sigma_params,
    est_gamma_params = est_gamma_params,
    est_nb2_params = est_nb2_params,
    gp_theta_prior = gp_theta_prior,
    use_expansion_prior = as.integer(expansion_prior)
  )

  if (is.null(par_list)) {
    pars <- c("x", "Z", "log_lik", "xstar")
    if (expansion_prior) pars <- c(pars, "psi")

    # family
    if (family %in% c("gaussian", "lognormal")) pars <- c(pars, "sigma")
    if (family %in% c("gamma")) pars <- c(pars, "gamma_a")
    if (family %in% c("nbinom2")) pars <- c(pars, "nb2_phi")

    if (est_correlation) pars <- c(pars, "Omega", "Sigma") # add correlation matrix
    if (estimate_nu) pars <- c(pars, "nu")
    if (estimate_trend_ar) pars <- c(pars, "phi")
    if (estimate_trend_ma) pars <- c(pars, "theta")
    if (!is.null(obs_covar)) pars <- c(pars, "b_obs")
    if (!is.null(pro_covar)) pars <- c(pars, "b_pro")
    if (est_sigma_process) pars <- c(pars, "sigma_process")
    if (trend_model == "gp") pars <- c(pars, "gp_theta")
    # if par list = "all", monitor everything --
    if (!is.null(par_list)) {
      if (par_list[1] == "all") {
        pars <- NA # removed pred
      }
    }
  } else {
    pars <- par_list
  }



  sampling_args <- list(
    object = stanmodels$dfa,
    data = data_list,
    pars = pars,
    control = control,
    chains = chains,
    iter = iter,
    thin = thin,
    show_messages = verbose,
    ...
  )

  out <- list()
  if (estimation[1]=="sampling") {
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
  }
  if(estimation[1]=="optimizing") {
    sampling_args <- list(
      object = stanmodels$dfa,
      data = data_list,
      verbose = verbose,
      ...
    )
    mod <- do.call(optimizing, sampling_args)
    out <- list(model = mod)
  }
  if(estimation[1]=="vb") {
    sampling_args <- list(
      object = stanmodels$dfa,
      data = data_list,
      iter = iter,
      pars = pars,
      ...
    )
    mod <- do.call(vb, sampling_args)
    out <- list(model = mod)
  }

  out[["sampling_args"]] <- sampling_args
  out[["orig_data"]] <- orig_data
  out[["shape"]] <- data_shape
  out[["z_model"]] <- z_model
  out[["z_bound"]] <- z_bound
  out[["trend_model"]] <- trend_model
  out[["estimation"]] <- estimation
  out[["scale"]] <- scale[1]
  out[["obs_covar"]] <- obs_covar
  out[["pro_covar"]] <- pro_covar
  out[["family"]] <- family

  out <- structure(out, class = "bayesdfa")
  out
}
