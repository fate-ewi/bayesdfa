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
#' @seealso invert_chains
#' @export
#'
#' @examples
#' \dontrun{
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' set.seed(1)
#' m <- fit_dfa(y = y, num_trends = 2, iter = 2000, chains = 4)
#' # chains were already inverted, but we can redo that, as an example, with:
#' find_inverted_chains(m$model, trend = 1, plot = TRUE)
#' }

find_inverted_chains <- function(model, trend = 1, thresh = 0.4, plot = FALSE) {
  e <- rstan::extract(model, permuted = FALSE)
  v <- reshape2::melt(e)

  vv <- v[grepl(paste0("x\\[", trend), v$parameters), ]
  # vv <- v[grepl(paste0("Z\\[[0-9]+,", trend, "\\]"), v$parameters), ]
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
    diff_ <- abs(out_df$mean_diff[1] - out_df$mean_diff[3])
    if (diff_ < thresh) {
      warning(paste0("Best permutation does not exceed threshold for trend ",
        trend, ". (", round(diff_, 2), " difference in mean SD)"))
    }
  }

  as.numeric(strsplit(out_df$neg[1], " ")[[1]])
}

#' Invert chains
#'
#' @param model A Stan model.
#' @param trends The number of trends in the DFA.
#' @param print Logical indicating whether the summary should be printed.
#' @param ... Other arguments to pass to \code{\link{find_inverted_chains}}.
#'
#' @export
#' @seealso find_inverted_chains
invert_chains <- function(model, trends = 1, print = FALSE, ...) {

  e <- rstan::extract(model, permuted = FALSE)
  ep <- rstan::extract(model, permuted = TRUE)
  pars <- colnames(e[1,,])
  n_mcmc <- dim(ep$Z)[1]
  n_chains <- dim(e)[2]
  ii <- c(seq(1, n_mcmc, n_mcmc/n_chains), n_mcmc + 1)

  for (k in seq_len(trends)) {
    f <- find_inverted_chains(model, trend = k)
    message(paste("Inverting chains", paste(f, collapse = " & "), "for trend", k))

    for (f_ in f) {
      for (i in grep(paste0("x\\[", k), pars)) {
        e[,f_,i] <- e[,f_,i] * -1
      }
      for (i in grep(paste0("Z\\[[0-9]+,", k, "\\]"), pars)) {
        e[,f_,i] <- e[,f_,i] * -1
      }
      # permuted
      # ep$Z[seq(ii[f_], ii[f_ +1] - 1),,k] <- ep$Z[seq(ii[f_], ii[f_ +1] - 1),,k] * -1
      # ep$x[seq(ii[f_], ii[f_ +1] - 1),,k] <- ep$x[seq(ii[f_], ii[f_ +1] - 1),,k] * -1
    }
  }
  # plot(ep$x[,,1])
  # plot(ep$x[,,2])

  mon <- rstan::monitor(e, print = print, warmup = 0)
  invisible(list(model = model, samples_permuted = ep, samples = e, monitor = mon))
}
