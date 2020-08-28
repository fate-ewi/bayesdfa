#' Plot the loadings from a DFA
#'
#' @param rotated_modelfit Output from [rotate_trends()].
#' @param names An optional vector of names for plotting the loadings.
#' @param facet Logical. Should there be a separate facet for each trend?
#'   Defaults to `TRUE`.
#' @param violin Logical. Should the full posterior densities be shown as a
#'   violin plot? Defaults to `TRUE`.
#' @param conf_level Confidence level for credible intervals. Defaults to 0.95.
#' @param threshold Numeric (0-1). Optional for plots, if included, only plot
#' loadings who have Pr(<0) or Pr(>0) > threshold. For example `threshold = 0.8`
#' would only display estimates where 80% of posterior density was above/below
#' zero. Defaults to `NULL` (not used).
#'
#' @seealso plot_trends fit_dfa rotate_trends
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point xlab ylab theme_bw theme aes_string
#'   element_blank position_dodge ggtitle geom_errorbar
#'   element_line element_text geom_line geom_violin coord_flip geom_hline
#'
#' @examples
#' set.seed(42)
#' s <- sim_dfa(num_trends = 2, num_ts = 4, num_years = 10)
#' # only 1 chain and 180 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, num_trends = 2, iter = 50, chains = 1)
#' r <- rotate_trends(m)
#' plot_loadings(r, violin = FALSE, facet = TRUE)
#' plot_loadings(r, violin = FALSE, facet = FALSE)
#' plot_loadings(r, violin = TRUE, facet = FALSE)
#' plot_loadings(r, violin = TRUE, facet = TRUE)

plot_loadings <- function(rotated_modelfit,
                          names = NULL,
                          facet = TRUE,
                          violin = TRUE,
                          conf_level = 0.95,
                          threshold=NULL) {
  v <- reshape2::melt(rotated_modelfit$Z_rot,
    varnames = c("iter", "name", "trend")
  )
  v$trend <- paste0("Trend ", v$trend)
  v$trend <- as.factor(v$trend)
  if (!is.null(names)) v$name <- names[v$name]
  v$name <- as.factor(v$name)

  v <- dplyr::group_by(v, .data$name, .data$trend)
  v <- dplyr::mutate(v,
    q_lower = sum(.data$value < 0) / length(.data$value),
    q_upper = 1 - .data$q_lower,
    prob_diff0 = max(.data$q_lower, .data$q_upper)
  )
  v <- dplyr::ungroup(v)

  vsum <- dplyr::group_by(v, .data$name, .data$trend)
  vsum <- dplyr::summarize(vsum,
    lower = quantile(.data$value, probs = (1 - conf_level) / 2),
    upper = quantile(.data$value, probs = 1 - (1 - conf_level) / 2),
    median = median(.data$value),
    q_lower = sum(.data$value < 0) / length(.data$value),
    q_upper = 1 - .data$q_lower,
    prob_diff0 = max(.data$q_lower, .data$q_upper)
  )
  df <- dplyr::ungroup(vsum)

  # filter values below threshold
  if (!is.null(threshold)) {
    df <- df[df$prob_diff0 >= threshold, ]
    v <- v[v$prob_diff0 >= threshold, ]
  }

  if (!violin) {
    p1 <- ggplot(df, aes_string(
      x = "name", y = "median", col = "trend",
      alpha = "prob_diff0"
    )) +
      geom_point(size = 3, position = position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
        position = position_dodge(0.3), width = 0
      ) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Time Series") + ylab("Loading")
  }

  if (violin) {
    p1 <- ggplot(v, aes_string(
      x = "name", y = "value", fill = "trend",
      alpha = "prob_diff0"
    )) +
      geom_violin(color = NA) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Time Series") + ylab("Loading")
  }

  if (facet) {
    p1 <- p1 + facet_wrap(~ trend, scales = "free_x")
  }

  p1
}
