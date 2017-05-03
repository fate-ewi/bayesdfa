#' Find outlying "black swan" jumps in trends
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}
#' @param threshold A probability threshold below which to
#'   flag trend events as extreme
#' @examples
#' \dontrun{
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays")])
#' y[,10] <- y[,10] + 2
#' set.seed(1)
#' m <- fit_dfa(y = y, num_trends = 2, iter = 1000, chains = 1, nu_fixed = 2)
#' r <- rotate_trends(m)
#' p <- plot_trends(r)
#' print(p)
#' find_swans(r)
#' }
#' @export

find_swans <- function(rotated_modelfit, threshold = 0.01) {
  x <- rotated_modelfit$trends_mean
  d <- apply(x, 1, function(xx) c(NA, diff(xx)))
  sds = apply(d, 2, sd, na.rm=T) # sds != 1

  prob = matrix(NA, nrow(d), ncol(d))
  for(i in 1:ncol(d)) {
    prob[,i] <- 1 - pnorm(abs(d[,i]), 0, sds[i])
  }
  #prob <- 1 - apply(d, 2, function(xx) pnorm(abs(xx), 0, 1))
  prob <- as.data.frame(prob)
  trends <- as.data.frame(t(x))
  trends$time <- seq_len(nrow(trends))
  prob$time <- seq_len(nrow(prob))
  trends <- reshape2::melt(trends, id.vars = c("time"))
  prob <- reshape2::melt(prob, id.vars = c("time"))
  names(trends) <- c("time", "trend_number", "trend_value")
  names(prob) <- c("time", "trend_number", "probability")

  trends <- dplyr::inner_join(trends, prob)
  trends$trend_number <- as.numeric(sub("V", "", trends$trend_number))
  trends$below_threshold <- trends$probability < threshold
  trends
  # ggplot(trends, aes(time, trend_value, color = below_threshold)) + geom_point() + facet_wrap(~trend_number)
  # ggplot(trends, aes(time, probability, color = below_threshold)) + geom_point() + facet_wrap(~trend_number)

}
# library(mvtnorm)

