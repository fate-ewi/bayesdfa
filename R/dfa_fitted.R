#' Get the fitted values from a DFA as a data frame
#'
#' @param modelfit Output from \code{\link{fit_dfa}}.
#' @param conf_level Probability level for CI.
#' @param names Optional vector of names for time series labels. Should be same length as the number of time series.
#'
#' @export
#' @return A data frame with the following columns: `ID` is an identifier for each time series, `time` is the time step, `y` is the observed values standardized to mean 0 and unit variance, `estimate` is the mean fitted value, `lower` is the lower CI, and `upper` is the upper CI.
#'
#' @seealso predicted plot_fitted fit_dfa
#'
#' @examples
#' \donttest{
#' y <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 4)
#' m <- fit_dfa(y = y$y_sim, num_trends = 2, iter = 50, chains = 1)
#' fitted <- dfa_fitted(m)
#' }
dfa_fitted <- function(modelfit, conf_level = 0.95, names = NULL) {

  # pred and Y have same dimensions if data is wide
  pred <- predicted(modelfit)
  n_mcmc <- dim(pred)[1]
  n_chains <- dim(pred)[2]
  n_years <- dim(pred)[3]
  n_ts <- dim(pred)[4]

  # this is the same for both data types
  df_pred <- data.frame(
    "ID" = rep(seq_len(n_ts), n_years),
    "time" = sort(rep(seq_len(n_years), n_ts)),
    "estimate" = c(t(apply(pred, c(3, 4), mean))),
    "lower" = c(t(apply(pred, c(3, 4), quantile, 1 - (1 - conf_level) / 2))),
    "upper" = c(t(apply(pred, c(3, 4), quantile, (1 - conf_level) / 2)))
  )

  if (modelfit$shape == "wide") {
    df_obs <- data.frame(
      "ID" = rep(seq_len(n_ts), n_years),
      "time" = sort(rep(seq_len(n_years), n_ts)),
      "y" = c(modelfit$orig_data)
    )
  } else {
    df_obs <- data.frame(
      "ID" = modelfit$orig_data[["ts"]],
      "time" = modelfit$orig_data[["time"]],
      "y" = modelfit$orig_data[["obs"]]
    )
  }
  df_obs$time <- df_obs$time - min(df_obs$time) + 1

  # standardize
  for (i in seq_len(n_ts)) {
    indx <- which(df_obs[["ID"]] == i)
    df_obs[indx, "y"] <- scale(df_obs[indx, "y"], center = TRUE, scale = TRUE)
  }

  df_obs <- df_obs[order(df_obs$ID, df_obs$time), ]
  df_pred <- df_pred[order(df_pred$ID, df_pred$time), ]

  if (!is.null(names)) {
    if (length(names) != n_ts) {
      warning("bayesdfa: Length of 'names' should match number of time series. Ignoring 'names'.")
    } else {
      df_pred$ID <- names[df_pred$ID]
      df_obs$ID <- names[df_obs$ID]
    }
  }

  df <- merge(df_obs, df_pred, by = c("ID", "time"), sort = FALSE)
  return(df)
}
