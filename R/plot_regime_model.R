#' Plot the state probabilities from [find_regimes()]
#'
#' @param model A model returned by [find_regimes()].
#' @param probs A numeric vector of quantiles to plot the credible intervals at.
#'   Defaults to `c(0.05, 0.95)`.
#' @param type Whether to plot the probabilities (default) or means.
#' @param regime_prob_threshold The probability density that must be above 0.5.
#'   Defaults to 0.9 before we classify a regime (only affects `"means"` plot).
#' @param plot_prob_indices Optional indices of probability plots to plot.
#'   Defaults to showing all.
#' @param flip_regimes Optional whether to flip regimes in plots, defaults to FALSE
#' @details Note that the original timeseries data (dots) are shown scaled
#'   between 0 and 1.
#' @importFrom dplyr "%>%"
#' @importFrom rlang .data
#' @export
#' @examples
#' \donttest{
#' data(Nile)
#' m <- fit_regimes(log(Nile), n_regimes = 2, chains = 1, iter = 50)
#' plot_regime_model(m)
#' plot_regime_model(m, plot_prob_indices=c(2))
#' plot_regime_model(m, type = "means")
#' }

plot_regime_model <- function(model, probs = c(0.05, 0.95),
                              type = c("probability", "means"),
                              regime_prob_threshold = 0.9,
                              plot_prob_indices = NULL,
                              flip_regimes = FALSE) {
  type <- match.arg(type)
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
    if (length(w) == 0) NA else w
  })

  if (flip_regimes) {
    mu_k <- 1 - mu_k
    u <- 1 - u
    l <- 1 - l
    med <- 1 - med
  }

  if (is.null(plot_prob_indices)) {
    # then plot all panels
    plot_prob_indices <- seq_len(ncol(med))
  }

  if (type == "probability") {
    df_l <- reshape2::melt(l, varnames = c("Time", "State"), value.name = "lwr")
    df_u <- reshape2::melt(u, varnames = c("Time", "State"), value.name = "upr")
    df_m <- reshape2::melt(med, varnames = c("Time", "State"), value.name = "median")
    df_y <- data.frame(y = range01(model$y), Time = seq_along(model$y))

    dplyr::inner_join(df_l, df_u, by = c("Time", "State")) %>%
      dplyr::inner_join(df_m, by = c("Time", "State")) %>%
      dplyr::filter(.data$State %in% plot_prob_indices) %>%
      dplyr::mutate(State = paste("State", .data$State)) %>%
      ggplot2::ggplot(
        ggplot2::aes_string("Time", y = "median", ymin = "lwr", ymax = "upr")
      ) +
      ggplot2::geom_ribbon(fill = "grey60") +
      ggplot2::geom_line(colour = "grey10", lwd = 0.8) +
      ggplot2::facet_wrap(~State) +
      ggplot2::coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
      ggplot2::geom_point(
        data = df_y,
        ggplot2::aes_string(x = "Time", y = "y"), inherit.aes = FALSE
      ) +
      ggplot2::ylab("Probability")
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
