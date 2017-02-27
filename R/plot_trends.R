#' Plot the trends from a DFA
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}
#' @param years Optional numeric vector of years for the plot
#'
#' @export
#'
#' @importFrom ggplot2 geom_ribbon facet_wrap
#'
#' @examples
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' m <- fit_dfa(y = y, num_trends = 2, iter = 1000)
#' r <- rotate_trends(m)
#' p <- plot_trends(r)
#' print(p)

plot_trends = function(rotated_modelfit, years=NULL) {
  # rotate the trends
  rotated = rotated_modelfit

  n_ts = dim(rotated$Z_rot)[2]
  n_trends = dim(rotated$Z_rot)[3]
  n_years = dim(rotated$trends_mean)[2]
  if(is.null(years)) years = 1:n_years

  # convert to df for ggplot
  df = data.frame(
    x = c(t(rotated$trends_mean)),
    lo = c(t(rotated$trends_lower)),
    hi = c(t(rotated$trends_upper)),
    trend = paste0("Trend ", sort(rep(1:n_trends, n_years))),
    time = rep(years, n_trends))

  # make faceted ribbon plot of trends
  p1 = ggplot(df, aes_string(x = "time", y = "x")) +
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) +
    geom_line() + facet_wrap("trend") +
    xlab("Time") + ylab("")
  p1

}
