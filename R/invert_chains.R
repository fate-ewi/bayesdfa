#' Find which chains to invert
#'
#' Find which chains to invert by checking the sum of the squared
#' deviations between the first chain and each other chain.
#'
#' @param model A Stan model, `rstanfit` object
#' @param trend Which trend to check
#' @param plot Logical: should a plot of the trend for each chain be made?
#'   Defaults to `FALSE`
#'
#' @importFrom ggplot2 geom_line
#' @importFrom utils combn
#' @seealso invert_chains
#' @export
#'
#' @examples
#' set.seed(2)
#' s <- sim_dfa(num_trends = 2)
#' set.seed(1)
#' m <- fit_dfa(y = s$y_sim, num_trends = 1, iter = 30, chains = 2)
#' # chains were already inverted, but we can redo that, as an example, with:
#' find_inverted_chains(m$model, plot = TRUE)

find_inverted_chains <- function(model, trend = 1, plot = FALSE) {

  chains = NULL # required for dplyr 0.8 update
  parameters = NULL
  value = NULL

  e <- rstan::extract(model, permuted = FALSE)
  v <- reshape2::melt(e)
  vv <- v[grepl(paste0("x\\[", trend), v$parameters), ]
  vv$parameters = as.factor(as.character(vv$parameters)) # needed with dplyr 0.8, all levels returned otherwise
  vv <- dplyr::group_by(vv, chains, parameters)
  vv <- dplyr::summarise(vv, estimate = stats::median(value))
  zz <- v[grepl(paste0("Z\\["), v$parameters), ]
  zz$parameters = as.factor(as.character(zz$parameters)) # needed with dplyr 0.8, all levels returned otherwise
  zz <- zz[grepl(paste0(trend, "]"), zz$parameters), ]
  zz <- dplyr::group_by(zz, chains, parameters)
  zz <- dplyr::summarise(zz, estimate = stats::median(value))
  ## vv is dimensioned nchains * nyears (x[1:nyears,trend=i])
  ## zz is dimensioned n_time series  (Z[1:time series,trend=i])

  if (plot) {
    p <- ggplot(vv, aes_string("as.numeric(parameters)", "estimate",
      color = "chains"
    )) +
      geom_line()
    print(p)
  }

  # cast parameters to df
  vvv <- reshape2::dcast(vv, parameters ~ chains, value.var = "estimate")
  vvv$parameters <- NULL
  zzz <- reshape2::dcast(zz, parameters ~ chains, value.var = "estimate")
  zzz$parameters <- NULL

  nchains <- ncol(vvv)

  # n_ts x n_years prediction matrix of product of trends and loadings
  flipped_chains <- 0
  pred0_loadings <- zzz[, 1] # loadings on first trend
  pred0_trend <- vvv[, 1] # loadings on second trend
  if (nchains > 1) {
    for (i in seq(2, nchains)) {
      pred1_loadings <- zzz[, i]
      pred1_trend <- vvv[, i]
      # see if flipped trend + loadings are more similar to chain 1 than not flipped
      if ((sum((-1 * pred1_loadings - pred0_loadings)^2) +
        sum((-1 * pred1_trend - pred0_trend)^2)) <
        (sum((pred1_loadings - pred0_loadings)^2) +
          sum((pred1_trend - pred0_trend)^2))) {
        # flip this chain -- seems to be something not right with commented out line
        #flipped_chains <- ifelse(flipped_chains == 0, i, c(flipped_chains, i))
        if(flipped_chains==0) {
          flipped_chains = i
        } else {
          flipped_chains = c(flipped_chains, i)
        }
      }
    }
  }
  flipped_chains
}

#' Invert chains
#'
#' @param model A Stan model, rstanfit object
#' @param trends The number of trends in the DFA, defaults to 1
#' @param print Logical indicating whether the summary should be printed.
#'   Defaults to `FALSE`.
#' @param ... Other arguments to pass to [find_inverted_chains()].
#'
#' @export
#' @seealso find_inverted_chains
invert_chains <- function(model, trends = 1, print = FALSE, ...) {
  e <- rstan::extract(model, permuted = FALSE)
  ep <- rstan::extract(model, permuted = TRUE)
  pars <- colnames(e[1, , ])
  n_mcmc <- dim(ep$Z)[1]
  n_chains <- dim(e)[2]

  for (k in seq_len(trends)) {
    f <- find_inverted_chains(model, trend = k)
    message(paste("Inverting chains", paste(f, collapse = " & "), "for trend", k))

    for (f_ in f) {
      for (i in grep(paste0("x\\[", k), pars)) {
        e[, f_, i] <- -1*e[, f_, i]
      }
      for (i in grep(paste0("Z\\[[0-9]+,", k, "\\]"), pars)) {
        e[, f_, i] <- -1*e[, f_, i]
      }
    }
  }

  mon <- rstan::monitor(e, print = print, warmup = 0)
  invisible(list(model = model, samples_permuted = ep, samples = e, monitor = mon))
}
