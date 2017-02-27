#' Plot the loadings from a DFA
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}.
#' @param names An optional vector of names for plotting the loadings.
#' @param threshold Numeric (0-1). Optional, if included, only plot loadings who have Pr(<0) or Pr(>0) > threshold
#' @param facet Logical. Should there be a separate facet for each trend?
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point xlab ylab theme_bw theme aes_string
#'   element_blank position_dodge ggtitle geom_errorbar
#'   element_line element_text geom_line
#'
#' @references
#' Del Negro, M., & Otrok, C. (2008). Dynamic factor models with time-varying
#' parameters: measuring changes in international business cycles.
#'
#' @examples
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])
#' m <- fit_dfa(y = y, num_trends = 2, iter = 500)
#' r <- rotate_trends(m)
#' p <- plot_loadings(r, threshold = 0.01)
#' print(p)

plot_loadings = function(rotated_modelfit,
  names = NULL,
  threshold = 0.8,
  facet = TRUE) {
  # rotate the trends
  rotated = rotated_modelfit

  n_ts = dim(rotated$Z_rot)[2]
  if (is.null(names)) {
    names = paste0("ts_", 1:n_ts)
  }
  n_trends = dim(rotated$Z_rot)[3]

  # convert to df for ggplot
  df = data.frame(
    x = c(rotated$Z_rot_mean),
    trend = paste0("Trend ", sort(rep(seq_len(n_trends), n_ts))),
    name = rep(names, n_trends),
    lower = c(rotated$Z_rot_mean) - c(apply(rotated$Z_rot,c(2,3),sd)),
    upper = c(rotated$Z_rot_mean) + c(apply(rotated$Z_rot,c(2,3),sd)),
    q_lower = c(apply(rotated$Z_rot,c(2,3),function(x) {return(length(which(x<0))/length(x))}))
  )
  df$q_upper = 1 - df$q_lower

  # replace low values with NAs
  df$x = ifelse((df$q_lower > threshold | df$q_upper > threshold), df$x, NA)
  df$lower = ifelse((df$q_lower > threshold | df$q_upper > threshold), df$lower, NA)
  df$upper = ifelse((df$q_lower > threshold | df$q_upper > threshold), df$upper, NA)

  # make faceted ribbon plot of trends
  if (facet) {
    #p1 = ggplot(df[!is.na(df$x),], aes_string(x = "name", y = "x")) +
    #  geom_point(position=position_dodge(0.3)) + facet_wrap("trend") +
    #  geom_errorbar(aes(ymin=lower, ymax=upper),alpha=0.5,position=position_dodge(0.3)) +
    #  xlab("Time Series") + ylab("Loading") +
    #  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #  theme(axis.text=element_text(size=5))

    loadings_plot = list()
    for(i in 1:n_trends) {
      loadings_plot[[i]] = ggplot(df[!is.na(df$x) & df$trend==paste0("Trend ",i),], aes_string(x = "name", y = "x")) +
        geom_point(position=position_dodge(0.3)) +
        geom_errorbar(aes_string(ymin="lower", ymax="upper"),alpha=0.5,position=position_dodge(0.3)) +
        xlab("") + ylab("Loading") + ggtitle(paste0("Trend ",i)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        theme(axis.text=element_text(size=5))
    }
    textstr = "gridExtra::grid.arrange(loadings_plot[[1]]"
    for(i in 2:n_trends) {
      # only show trends with > 1 loading
      if(nrow(df[!is.na(df$x) & df$trend==paste0("Trend ",i),]) > 0) {
      textstr = paste0(textstr, ", loadings_plot[[",i,"]]")
      }
    }
    textstr = paste0(textstr, ")")
    p1 = eval(parse(text=textstr))
  }
  if (!facet) {
    p1 = ggplot(df[!is.na(df$x),], aes_string(x = "name", y = "x", col = "trend")) +
      geom_point(size = 3, alpha = 0.5,position=position_dodge(0.3)) +
      geom_errorbar(aes_string(ymin="lower", ymax="upper"),alpha=0.5,position=position_dodge(0.3)) +
      xlab("Time Series") + ylab("Loading") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(axis.text=element_text(size=5))
  }
  p1
}
