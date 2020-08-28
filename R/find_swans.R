#' Find outlying "black swan" jumps in trends
#'
#' @param rotated_modelfit Output from [rotate_trends()].
#' @param threshold A probability threshold below which to flag trend events as
#'   extreme
#' @param plot Logical: should a plot be made?
#'
#' @return
#' Prints a ggplot2 plot if `plot = TRUE`; returns a data frame indicating the
#' probability that any given point in time represents a "black swan" event
#' invisibly.
#'
#' @examples
#' set.seed(1)
#' s <- sim_dfa(num_trends = 1, num_ts = 3, num_years = 30)
#' s$y_sim[1, 15] <- s$y_sim[1, 15] - 6
#' plot(s$y_sim[1,], type = "o")
#' abline(v = 15, col = "red")
#' # only 1 chain and 250 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, num_trends = 1, iter = 50, chains = 1, nu_fixed = 2)
#' r <- rotate_trends(m)
#' p <- plot_trends(r) #+ geom_vline(xintercept = 15, colour = "red")
#' print(p)
#' # a 1 in 1000 probability if was from a normal distribution:
#' find_swans(r, plot = TRUE, threshold = 0.001)
#'
#' @references
#' Anderson, S.C., Branch, T.A., Cooper, A.B., and Dulvy, N.K. 2017.
#' Black-swan events in animal populations. Proceedings of the National Academy
#' of Sciences 114(12): 3252–3257. https://doi.org/10.1073/pnas.1611525114
#'
#' @export
#' @importFrom stats pnorm

find_swans <- function(rotated_modelfit,
  threshold = 0.01,
  plot = FALSE) {

  x <- rotated_modelfit$trends_mean
  d <- apply(x, 1, function(xx) c(NA, diff(xx)))
  sds <- apply(d, 2, sd, na.rm = TRUE) # sds != 1

  prob <- matrix(NA, nrow(d), ncol(d))
  for (i in seq_len(ncol(d))) {
    prob[, i] <- 1 - pnorm(abs(d[, i]), 0, sds[i])
  }
  prob <- as.data.frame(prob)
  trends <- as.data.frame(t(x))
  trends$time <- seq_len(nrow(trends))
  prob$time <- seq_len(nrow(prob))
  trends <- reshape2::melt(trends, id.vars = c("time"))
  prob <- reshape2::melt(prob, id.vars = c("time"))
  names(trends) <- c("time", "trend_number", "trend_value")
  names(prob) <- c("time", "trend_number", "probability")

  trends$trend_number <- as.character(sub("V", "", trends$trend_number))
  prob$trend_number <- as.character(sub("V", "", prob$trend_number))
  trends <- dplyr::inner_join(trends, prob, c("time", "trend_number"))
  trends$below_threshold <- trends$probability < threshold

  if (plot) {
    g <- ggplot(trends, aes_string(
      x = "time", y = "trend_value",
      color = "below_threshold"
    )) +
      geom_point() + facet_wrap(~ trend_number)
    print(g)
  }
  invisible(trends)
}
