#' Plot the loadings from a DFA
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}.
#' @param names An optional vector of names for plotting the loadings.
#' @param threshold Threshold above which to plot the loadings.
#' @param facet Logical. Should there be a separate facet for each trend?
#'
#' @export
#'
#' @importFrom ggplot2 ggplot geom_point xlab ylab theme_bw theme aes_string
#'   element_blank
#'   element_line element_text geom_line

plot_loadings = function(rotated_modelfit, names=NULL, threshold=0.1, facet=TRUE) {
  # rotate the trends
  rotated = rotated_modelfit

  n_ts = dim(rotated$Z_rot)[2]
  if(is.null(names)) {
    names=paste0("ts_",1:n_ts)
  }
  n_trends = dim(rotated$Z_rot)[3]

  # convert to df for ggplot
  df = data.frame("x" = c(rotated$Z_rot_mean),
    "trend" = paste0("Trend ",sort(rep(1:n_trends,n_ts))),
    "name"=rep(names,n_trends))

  # replace low values with NAs
  df$x = ifelse(abs(df$x) < threshold, NA, df$x)

  # make faceted ribbon plot of trends
  if(facet==TRUE) {
  p1 = ggplot(df, aes_string(x = "name", y = "x")) +
    geom_point() + facet_wrap("trend") +
    xlab("Time Series") + ylab("Loading") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }
  if(facet==FALSE) {
    p1 = ggplot(df, aes_string(x = "name", y = "x", col="trend")) +
      geom_point(size=3,alpha=0.5) +
      xlab("Time Series") + ylab("Loading") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }
  print(p1)
}
