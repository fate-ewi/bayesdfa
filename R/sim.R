#' Simulate from a DFA
#'
#' @param num_trends The number of trends
#' @param num_years The number of years
#' @param num_ts The number of timeseries
#' @param loadings_matrix A loadings matrix. The number of rows should match the
#'   number of timeseries and the number of columns should match the number of
#'   trends. Note that this loadings matrix will be internally manipulated by
#'   setting some elements to 0 and constraining some elements to 1 so that the
#'   model can be fitted. See \code{\link{fit_dfa}}. See the outfit element
#'   \code{Z} in the returned list is to see the manipulated loadings matrix.
#' @param sigma A vector of standard deviations on the observation error. Should
#'   be of the same length as the number of trends.
#'
#' @return A list with the following elements: y_sim is the simulated data, pred
#'   is the true underlying data without observation error added, x is the
#'   underlying trends, Z is the manipulated loadings matrix that is fed to the
#'   model.

sim_dfa <- function(
  num_trends = 2,
  num_years = 20,
  num_ts = 4,
  loadings_matrix = matrix(nrow = num_ts, ncol = num_trends,
    rnorm(num_ts * num_trends, 0, 1)),
  sigma = rlnorm(num_trends, meanlog = log(0.2), 0.1)
) {

  y_ignore <- matrix(rnorm(num_ts * num_years), nrow = num_ts, ncol = num_years)
  d <- fit_dfa(y_ignore, num_trends = num_trends, sample = FALSE, zscore = FALSE)

  Z <- matrix(nrow = d$P, ncol = d$K)
  y <- vector(mode = "numeric", length = d$N)
  z <- as.numeric(loadings_matrix)

  for(i in seq_len(d$nZ)) {
    Z[d$row_indx[i],d$col_indx[i]] <- z[i];
  }

  # fill in zero elements
  if(d$nZero > 2) {
    for(i in seq_len(d$nZero-2)) {
      Z[d$row_indx_z[i],d$col_indx_z[i]] <- 0;
    }
  }
  for (k in seq_len(d$K)) {
    Z[k, k] <- 1 # add constraint for Z diagonal
  }

  x <- matrix(nrow = d$K, ncol = d$N) # random walk-trends

  # initial state for each trend
  for (k in seq_len(d$K)) {
    x[k, 1] <- rnorm(1, 0, 1)
    for (t in 2:d$N) {
      x[k, t] <- x[k, t - 1] + rt(1, df = d$nu_fixed) # random walk
    }
  }
  pred <- Z %*% x
  for (i in seq_len(d$n_pos)) {
    y[i] <- rnorm(1, pred[d$row_indx_pos[i], d$col_indx_pos[i]],
      sigma[d$varIndx[d$row_indx_pos[i]]])
  }
  y_sim <- matrix(y, nrow = d$P)
  list(y_sim = y_sim, pred = pred, x = x, Z = Z)
}
