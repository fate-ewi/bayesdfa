#' Plot the loadings from a DFA
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}.
#' @param names An optional vector of names for plotting the loadings.
#' @param threshold Numeric (0-1). Optional, if included, only plot loadings who
#'   have Pr(<0) or Pr(>0) > threshold
#' @param facet Logical. Should there be a separate facet for each trend?
#' @param violin Logical. Should the full posterior densities be shown as a
#'   violin plot?
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point xlab ylab theme_bw theme aes_string
#'   element_blank position_dodge ggtitle geom_errorbar
#'   element_line element_text geom_line geom_violin coord_flip geom_hline
#'
#' @references
#' Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
#' parameters: measuring changes in international business cycles.
#'
#' @examples
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' m <- fit_dfa(y = y, num_trends = 2, iter = 2000)
#' r <- rotate_trends(m)
#' plot_loadings(r, threshold = 0.01, violin = FALSE, facet = TRUE)
#' plot_loadings(r, threshold = 0.01, violin = FALSE, facet = FALSE)
#' plot_loadings(r, threshold = 0.01, violin = TRUE, facet = FALSE)
#' plot_loadings(r, threshold = 0.01, violin = TRUE, facet = TRUE)

plot_loadings = function(rotated_modelfit,
  names = NULL,
  threshold = 0.8,
  facet = TRUE,
  violin = TRUE,
  conf_level = 0.95) {

  v <- reshape2::melt(rotated_modelfit$Z_rot, varnames = c("iter", "name", "trend"))
  v$trend <- paste0("Trend ", v$trend)
  v$trend <- as.factor(v$trend)
  if (!is.null(names)) v$name <- names[v$name]
  v$name <- as.factor(v$name)

  vsum <- dplyr::group_by_(v, "name", "trend")
  vsum <- dplyr::summarize(vsum,
    lower = quantile(value, probs = (1 - conf_level) / 2),
    upper = quantile(value, probs = 1 - (1 - conf_level) / 2),
    median = median(value))
  df <- dplyr::ungroup(vsum)

  # q_lower = c(apply(rotated$Z_rot,c(2,3),function(x) {return(length(which(x<0))/length(x))}))
  # df$q_upper <- 1 - df$q_lower

  # replace low values with NAs
  # df$x <- ifelse((df$q_lower > threshold | df$q_upper > threshold), df$x, NA)
  # df$lower <- ifelse((df$q_lower > threshold | df$q_upper > threshold), df$lower, NA)
  # df$upper <- ifelse((df$q_lower > threshold | df$q_upper > threshold), df$upper, NA)

  # df <- df[!is.na(df$x), ]

  if (!violin) {
    p1 <- ggplot(df, aes_string(x = "name", y = "median", col = "trend")) +
      geom_point(size = 3, alpha = 0.95, position = position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin = "lower", ymax = "upper"),
        alpha = 0.6, position = position_dodge(0.3), width = 0) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Time Series") + ylab("Loading")
  }

  if (violin) {
    p1 <- ggplot(v, aes_string(x = "name", y = "value", fill = "trend")) +
      geom_violin(color = NA) +
      geom_hline(yintercept = 0, lty = 2) +
      coord_flip() + xlab("Time Series") + ylab("Loading")
  }

  if (facet)
    p1 <- p1 + facet_wrap(~trend, scales = "free_x")

  p1
}
