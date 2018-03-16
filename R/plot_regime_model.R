#' Plot the state probabilities from [find_regimes()]
#'
#' @param model A model returned by [find_regimes()].
#' @param probs A numeric vector of quantiles to plot the credible intervals at.
#'   Defaults to `c(0.05, 0.95)`.
#' @param type Whether to plot the probabilities (default) or means.
#' @param regime_prob_threshold The probability density that must be above 0.5.
#'   Defaults to 0.9 before we classify a regime (only affects `"means"` plot).
#' @details Note that the original timeseries data (dots) are shown scaled
#'   between 0 and 1.
#' @export
#' @examples
#' data(Nile)
#' m <- fit_regimes(log(Nile), n_regimes = 2, chains = 1, iter = 800)
#' plot_regime_model(m)
#' plot_regime_model(m, type = "means")

plot_regime_model <- function(model, probs = c(0.05, 0.95),
                              type = c("probability", "means"),
                              regime_prob_threshold = 0.9) {
  gamma_tk <- rstan::extract(model$model, pars = "gamma_tk")[[1]]
  mu_k <- rstan::extract(model$model, pars = "mu_k")[[1]]
  l <- apply(gamma_tk, 2:3, quantile, probs = probs[[1]])
  u <- apply(gamma_tk, 2:3, quantile, probs = probs[[2]])
  med <- apply(gamma_tk, 2:3, quantile, probs = 0.5)
  range01 <- function(x) (x - min(x)) / (max(x) - min(x))
  mu_k_low <- apply(mu_k, 2, quantile, probs = probs[[1]])
  mu_k_high <- apply(mu_k, 2, quantile, probs = probs[[2]])
  mu_k <- apply(mu_k, 2, median)
  confident_regimes <- apply(gamma_tk, 2:3, function(x)
    mean(x > 0.5) > regime_prob_threshold)
  regime_indexes <- apply(confident_regimes, 1, function(x) {
    w <- which(x)
    ifelse(length(w) == 0, NA, w)
  })

  if (type[[1]] == "probability") {
    oldpar <- par("mfrow")
    par(mfrow = c(1, ncol(med)))
    for (i in seq_len(ncol(med))) {
      plot(l[, i],
        ylim = c(0, 1), col = "grey40", lty = 2, type = "n",
        main = paste("State", LETTERS[i]), ylab = "Probability of being in given state",
        xlab = "Time"
      )
      polygon(c(1:nrow(u), nrow(u):1), c(l[, i], rev(u[, i])), col = "grey70", border = "grey70")
      lines(1:nrow(u), med[, i], col = "black", lwd = 2)
      points(1:nrow(u), range01(model$y), col = "#FF000070", pch = 3)
    }
    par(mfrow = oldpar)
  } else {
    plot(as.numeric(model$y),
      col = "#FF000070", pch = 3, ylab = "Time series value",
      xlab = "Time"
    )
    if (!all(is.na(regime_indexes))) {
      for (i in seq_along(regime_indexes)) {
        segments(
          x0 = i - 0.5, x1 = i + 0.5, y0 = mu_k[regime_indexes[i]],
          y1 = mu_k[regime_indexes[i]]
        )
        polygon(c(i - 0.5, i - 0.5, i + 0.5, i + 0.5),
          c(
            mu_k_low[regime_indexes[i]], mu_k_high[regime_indexes[i]],
            mu_k_high[regime_indexes[i]], mu_k_low[regime_indexes[i]]
          ),
          border = NA, col = "#00000050"
        )
      }
    }
  }
}
