#' Plot the fitted values from a DFA
#'
#' @param modelfit Output from \code{\link{fit_dfa}}, a rstanfit object
#' @param conf_level Probability level for CI.
#' @param names Optional vector of names for plotting labels TODO. Should be same length as the number of time series
#' @param spaghetti Defaults to FALSE, but if TRUE puts all raw time series (grey) and fitted values on a single plot
#'
#' @export
#' @seealso plot_loadings fit_dfa rotate_trends dfa_fitted
#'
#' @importFrom ggplot2 geom_ribbon facet_wrap scale_color_manual
#' @importFrom viridisLite viridis
#'
#' @examples
#' \donttest{
#' y <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 4)
#' m <- fit_dfa(y = y$y_sim, num_trends = 2, iter = 50, chains = 1)
#' p <- plot_fitted(m)
#' print(p)
#'
#' p <- plot_fitted(m, spaghetti = TRUE)
#' print(p)
#' }
plot_fitted <- function(modelfit, conf_level = 0.95, names = NULL, spaghetti = FALSE) {
  df <- dfa_fitted(modelfit, conf_level = conf_level, names = names)
  df$ID <- as.factor(df$ID)

  if (spaghetti == TRUE) {
    cols <- viridis(length(unique((df$ID))), end = 0.8)
    p1 <- ggplot(df) +
      geom_line(aes_string(x = "time", y = "y", group = "ID"),
        color = "grey50", size = 0.5
      ) +
      geom_line(aes_string(x = "time", y = "estimate", group = "ID", color = "ID"),
        size = 1.2
      ) +
      scale_color_manual(values = cols) +
      xlab("Time") +
      theme(legend.position = "none")
  } else {
    p1 <- ggplot(df) +
      geom_ribbon(aes_string(x = "time", ymin = "lower", ymax = "upper"), alpha = 0.4) +
      geom_line(aes_string(x = "time", y = "estimate")) +
      geom_point(aes_string(x = "time", y = "y"),
        col = "red",
        size = 0.5,
        alpha = 0.4
      ) +
      facet_wrap("ID", scales = "free_y") +
      xlab("Time") +
      ylab("")
  }
  p1
}
