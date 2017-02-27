#' Find which chains to invert
#'
#' Finds which chains to invert by finding the best combination of chains to
#' multiply by -1 to minimize the standard deviation of the trends across the
#' chains.
#'
#' @param model A Stan model
#' @param trend Which trend to check
#' @param thresh A threshold for the required absolute difference in mean
#'   standard deviation between the best permutation and the next best
#'   permutation.
#' @param plot Should a plot of the trend for each chain be made?
#'
#' @importFrom ggplot2 geom_line
#' @importFrom utils combn
#' @export
#'
#' @examples
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' set.seed(1)
#' m <- fit_dfa(y = y, num_trends = 3, iter = 600, chains = 4)
#' find_flipped_chains(m, trend = 1, plot = TRUE)

find_flipped_chains <- function(model, trend = 1, thresh = 0.8, plot = FALSE) {
  e <- extract(model, permuted = FALSE)
  v <- reshape2::melt(e)

  vv <- v[grepl(paste0("x\\[", trend), v$parameters), ]
  vv <- dplyr::group_by_(vv, "chains", "parameters")
  vv <- dplyr::summarize_(vv, estimate = "stats::median(value)")

  if (plot) {
    p <- ggplot(vv, aes_string("as.numeric(parameters)", "estimate",
      color = "chains")) +
      geom_line()
    print(p)
  }

  vvv <- reshape2::dcast(vv, parameters ~ chains, value.var = "estimate")
  vvv$parameters <- NULL

  nchains <- ncol(vvv)

  check <- combn(seq_len(nchains), 1, simplify = FALSE)

  if (nchains > 1) {
    for (i in seq(2, nchains)) {
      check <- c(check, combn(seq_len(nchains), i, simplify = FALSE))
    }
  }

  mean_diff <- vector(mode = "numeric", length = length(check))

  out_df <- data.frame(mean_diff)
  out_df$neg <- NA

  for (k in seq_along(check)) {
    rs <- vector(mode = "numeric", length = nchains^2)
    l <- 1
    x_temp <- vvv
    for (q in seq_along(check[[k]])) {
      x_temp[check[[k]][q]] <- -1 * x_temp[check[[k]][q]]
    }
    out_df[k, "mean_diff"] <- mean(apply(x_temp, 1, sd))
    out_df[k, "neg"] <- paste(check[[k]], collapse = " ")
  }

  out_df$nchar <- nchar(out_df$neg)
  out_df <- dplyr::arrange_(out_df, "mean_diff", "nchar")

  if (nrow(out_df) > 1) {
    if (abs(out_df$mean_diff[1] - out_df$mean_diff[3]) < thresh) {
      warning("Best permutation does not exceed threshold!")
    }
  }

  as.numeric(strsplit(out_df$neg[1], " ")[[1]])
}
#
# flip_chains <- function(model, chains_to_invert = c(1)) {
#
# }
