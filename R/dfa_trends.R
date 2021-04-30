#' Get the trends from a DFA as a data frame
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}.
#' @param years Optional numeric vector of years.
#'
#' @export
#' @return A data frame with the following columns: `time` is the time step, `trend_number` is an identifier for each trend, `estimate` is the trend mean, `lower` is the lower CI, and `upper` is the upper CI.
#'
#' @seealso plot_trends fit_dfa rotate_trends
#'
#' @examples
#' set.seed(1)
#' s <- sim_dfa(num_trends = 1)
#' m <- fit_dfa(y = s$y_sim, num_trends = 1, iter = 50, chains = 1)
#' r <- rotate_trends(m)
#' trends <- dfa_trends(r)
dfa_trends <- function(rotated_modelfit, years = NULL) {

  rotated <- rotated_modelfit
  n_ts <- dim(rotated$Z_rot)[2]
  n_trends <- dim(rotated$Z_rot)[3]

  n_years <- dim(rotated$trends_mean)[2]
  if (is.null(years)) years <- seq_len(n_years)

  df <- data.frame(
    time = rep(years, n_trends),
    trend_number = paste0("Trend ", sort(rep(seq_len(n_trends), n_years))),
    estimate = c(t(rotated$trends_mean)),
    lower = c(t(rotated$trends_lower)),
    upper = c(t(rotated$trends_upper)))

  return(df)
}
