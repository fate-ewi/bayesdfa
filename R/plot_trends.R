#' Plot the trends from a DFA
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}
#' @param years Optional numeric vector of years for the plot
#' @param highlight_outliers Logical. Should trend events
#'   that exceed the probability of occurring with a normal distribution as
#'   defined by \code{threshold} be highlighted? Defaults to FALSE
#' @param threshold A probability threshold below which to
#'   flag trend events as extreme. Defaults to 0.01
#'
#' @export
#' @seealso plot_loadings fit_dfa rotate_trends
#'
#' @importFrom ggplot2 geom_ribbon facet_wrap geom_point
#'
#' @examples
#' set.seed(1)
#' s <- sim_dfa(num_trends = 1)
#' m <- fit_dfa(y = s$y_sim, num_trends = 1, iter = 50, chains = 1)
#' r <- rotate_trends(m)
#' p <- plot_trends(r)
#' print(p)

plot_trends <- function(rotated_modelfit,
  years = NULL,
  highlight_outliers = FALSE,
  threshold = 0.01) {

  # rotate the trends
  rotated <- rotated_modelfit

  n_ts <- dim(rotated$Z_rot)[2]
  n_trends <- dim(rotated$Z_rot)[3]

  n_years <- dim(rotated$trends_mean)[2]
  if (is.null(years)) years <- seq_len(n_years)

  # convert to df for ggplot
  df <- data.frame(
    x = c(t(rotated$trends_mean)),
    lo = c(t(rotated$trends_lower)),
    hi = c(t(rotated$trends_upper)),
    trend = paste0("Trend ", sort(rep(seq_len(n_trends), n_years))),
    time = rep(years, n_trends)
  )

  # make faceted ribbon plot of trends
  p1 <- ggplot(df, aes_string(x = "time", y = "x")) +
    geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) +
    geom_line() + facet_wrap("trend") +
    xlab("Time") + ylab("")

  if (highlight_outliers) {
    swans <- find_swans(rotated, threshold = threshold)
    df$outliers <- swans$below_threshold
    p1 <- p1 + geom_point(data = df[which(df$outliers), ], color = "red")
  }

  p1
}
