#' Plot the trends from a DFA
#'
#' @param modelfit Output from \code{\link{fit_dfa}}, a rstanfit object
#' @param names Optional vector of names for plotting labels
#' @param spaghetti Defaults to FALSE, but if TRUE puts all raw time series (grey) and trends on a single plot
#'
#' @export
#' @seealso plot_loadings fit_dfa rotate_trends
#'
#' @importFrom ggplot2 geom_ribbon facet_wrap scale_color_identity
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

plot_fitted <- function(modelfit, names = NULL, spaghetti = FALSE) {
  if(modelfit$shape == "wide") {
    n_ts <- dim(modelfit$orig_data)[1]
    n_years <- dim(modelfit$orig_data)[2]
  } else {
    n_ts = max(modelfit$orig_data[,"ts"])
    n_years = max(modelfit$orig_data[,"time"])
  }

  # pred and Y have same dimensions if data is wide
  pred <- predicted(modelfit)

  # this is the same for both data types
  df_pred <- data.frame(
    "ID" = rep(seq_len(n_ts), n_years),
    "Time" = sort(rep(seq_len(n_years), n_ts)),
    "mean" = c(t(apply(pred, c(3, 4), mean))),
    "lo" = c(t(apply(pred, c(3, 4), quantile, 0.025))),
    "hi" = c(t(apply(pred, c(3, 4), quantile, 0.975)))
  )
  if(modelfit$shape == "wide") {
    df_obs <- data.frame(
    "ID" = rep(seq_len(n_ts), n_years),
    "Time" = sort(rep(seq_len(n_years), n_ts)),
    "y" = c(modelfit$orig_data))
  } else {
    df_obs <- data.frame(
      "ID" = modelfit$orig_data[,"ts"],
      "Time" = modelfit$orig_data[,"time"],
      "y" = modelfit$orig_data[,"obs"])
  }
  # standardize
  for(i in seq_len(n_ts)) {
    indx = which(df_obs[["ID"]] == i)
    df_obs[indx,"y"] = scale(df_obs[indx,"y" ], center = TRUE, scale = TRUE)
  }

  if (!is.null(names)) {
    df$ID <- names[df$ID]
  }

  if(spaghetti==TRUE) {

    cols = viridis(max(df_pred$ID), end=0.8)
    #hues = seq(15, 375, length = max(df_pred$ID) + 1)
    #cols = hcl(h = hues, l = 65, c = 100)[1:max(df_pred$ID)]
    df_obs$col = "grey50"
    df_obs$size=0.5
    names(df_pred)[which(names(df_pred)=="mean")]="y"
    df_pred$col = cols[df_pred$ID]
    df_pred$size=1.2
    df_pred$ID = df_pred$ID + max(df_obs$ID)

    df = rbind(df_obs,df_pred[,c("ID","Time","y","col","size")])
    df$ID = as.factor(df$ID)

    p1 <- ggplot(data=df, aes_string(x = "Time", y = "y", group="ID",color="col")) +
      geom_line(size=df$size) +
      scale_color_identity()
  } else {
    # make faceted ribbon plot of trends
    p1 <- ggplot(df_pred, aes_string(x = "Time", y = "mean")) +
      geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) +
      geom_line() +
      geom_point(data = df_obs,
        aes_string(x = "Time", y = "y"),
        col = "red",
        size = 0.5,
        alpha = 0.4
      ) +
      facet_wrap("ID", scales = "free_y") +
      xlab("Time") + ylab("")
  }
  p1
}
