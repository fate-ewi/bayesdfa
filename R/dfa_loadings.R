#' Get the loadings from a DFA as a data frame
#'
#' @param rotated_modelfit Output from \code{\link{rotate_trends}}.
#' @param names An optional vector of names for plotting the loadings.
#' @param summary Logical. Should the full posterior densities be returned? Defaults to `TRUE`.
#' @param conf_level Confidence level for credible intervals. Defaults to 0.95.
#'
#' @seealso plot_loadings fit_dfa rotate_trends
#'
#' @export
#' @return A data frame with the following columns:
#' `name` is an identifier for each loading, `trend` is the trend for the
#' loading, `median` is the posterior median loading, `lower` is the lower CI,
#' `upper` is the upper CI, and `prob_diff0` is the probability the loading is
#' different than 0. When `summary = FALSE`, there is no `lower` or `upper`
#' columns and instead there are columns `chain` and `draw`.
#'
#' @examples
#' set.seed(42)
#' s <- sim_dfa(num_trends = 2, num_ts = 4, num_years = 10)
#' # only 1 chain and 180 iterations used so example runs quickly:
#' m <- fit_dfa(y = s$y_sim, num_trends = 2, iter = 50, chains = 1)
#' r <- rotate_trends(m)
#' loadings <- dfa_loadings(r, summary = TRUE)
#' loadings <- dfa_loadings(r, summary = FALSE)
dfa_loadings <- function(rotated_modelfit,
                         names = NULL,
                         summary = TRUE,
                         conf_level = 0.95) {
  v <- reshape2::melt(rotated_modelfit$Z_rot,
    varnames = c("iter", "name", "trend"), value.name = "loading"
  )
  v$draw <- as.numeric(gsub("_chain.*$", "", v$iter))
  v$chain <- as.numeric(gsub("^[0-9]+_chain:", "", v$iter))
  v$iter <- NULL
  v <- v[, c("chain", "draw", "name", "trend", "loading")]

  v$trend <- paste0("Trend ", v$trend)
  v$trend <- as.factor(v$trend)
  if (!is.null(names)) v$name <- names[v$name]
  v$name <- as.factor(v$name)

  ## q_lower = proportion of draws less than zero
  ## q_upper = proportion of draws greater than zero
  v <- dplyr::group_by(v, .data$name, .data$trend)
  v <- dplyr::mutate(v,
    q_lower = sum(.data$loading < 0) / length(.data$loading),
    q_upper = 1 - .data$q_lower,
    prob_diff0 = max(.data$q_lower, .data$q_upper)
  )
  v <- as.data.frame(dplyr::ungroup(v))
  out <- v

  if (summary) {
    vsum <- dplyr::group_by(v, .data$name, .data$trend)
    vsum <- dplyr::summarize(vsum,
      median = median(.data$loading),
      lower = quantile(.data$loading, probs = (1 - conf_level) / 2),
      upper = quantile(.data$loading, probs = 1 - (1 - conf_level) / 2),
      q_lower = sum(.data$loading < 0) / length(.data$loading),
      q_upper = 1 - .data$q_lower,
      prob_diff0 = max(.data$q_lower, .data$q_upper)
    )
    df <- as.data.frame(dplyr::ungroup(vsum))
    out <- df
  }

  out$q_lower <- NULL
  out$q_upper <- NULL

  return(out)
}
