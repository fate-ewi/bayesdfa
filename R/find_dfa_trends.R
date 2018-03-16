#' Find the best number of trends according to LOOIC
#'
#' Fit a DFA with different number of trends and return the leave one out (LOO)
#' value as calculated by the [loo][loo::loo-package()] package.
#'
#' @param y A matrix of data to fit. Columns represent time element.
#' @param kmin Minimum number of trends.
#' @param kmax Maximum number of trans.
#' @param iter Iterations when sampling from each Stan model.
#' @param compare_normal If `TRUE`, does model selection comparison of Normal vs.
#'   Student-t errors
#' @param convergence_threshold The maximum allowed value of Rhat to determine
#'   convergence of parameters
#' @param variance Vector of variance arguments for searching over large groups
#'   of models. Can be either or both of ("equal","unequal")
#' @param ... Other arguments to pass to `fit_dfa()`
#' @export
#' @examples
#' \dontrun{
#' y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays")])
#' set.seed(1)
#' m <- find_dfa_trends(y = y, kmin = 1, kmax = 2,
#'   iter = 1000, chains = 1)
#' }
#' @importFrom loo loo extract_log_lik
#' @importFrom stats quantile time varimax
#' @importFrom rlang .data

find_dfa_trends <- function(y = y, kmin = 1, kmax = 5, iter = 2000,
                            compare_normal = FALSE, convergence_threshold = 1.05,
                            variance = c("equal", "unequal"), ...) {
  df <- data.frame(
    model = seq(1, ifelse(compare_normal == FALSE,
      length(variance) * length(seq(kmin, kmax)),
      2 * length(variance) * length(seq(kmin, kmax))
    )),
    num_trends = NA,
    looic = NA,
    cor = NA,
    error = NA,
    converge = FALSE,
    stringsAsFactors = FALSE
  )
  best_model <- NULL
  best_loo <- 1.0e50

  indx <- 1

  if (length(which(variance %in% "equal")) > 0) {
    for (i in seq(kmin, kmax)) {
      model <- fit_dfa(
        y = y, num_trends = i, iter = iter,
        estimate_nu = TRUE, ...
      )

      df$converge[indx] <- is_converged(model, convergence_threshold)
      df$num_trends[indx] <- i
      df$looic[indx] <- loo::loo(loo::extract_log_lik(model$model))$looic

      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
        best_model <- model
        best_loo <- df$looic[indx]
      }
      df$error[indx] <- "student-t"
      df$cor[indx] <- "equal"
      indx <- indx + 1
    }
  }

  if (length(which(variance %in% "unequal")) > 0) {
    for (i in seq(kmin, kmax)) {
      model <- fit_dfa(
        y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)),
        estimate_nu = TRUE, ...
      )
      df$num_trends[indx] <- i

      df$looic[indx] <- loo::loo(loo::extract_log_lik(model$model))$looic
      df$converge[indx] <- is_converged(model, convergence_threshold)
      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
        best_model <- model
        best_loo <- df$looic[indx]
      }
      df$error[indx] <- "student-t"
      df$cor[indx] <- "independent"
      indx <- indx + 1
    }
  }


  if (compare_normal == TRUE) {
    if (length(which(variance %in% "equal")) > 0) {
      for (i in seq(kmin, kmax)) {
        model <- fit_dfa(
          y = y, num_trends = i, iter = iter, nu_fixed = 100,
          estimate_nu = FALSE, ...
        )
        df$num_trends[indx] <- i
        df$looic[indx] <- loo::loo(loo::extract_log_lik(model$model))$looic

        df$converge[indx] <- is_converged(model, convergence_threshold)
        # if model is best, keep it
        if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
          best_model <- model
          best_loo <- df$looic[indx]
        }
        df$error[indx] <- "normal"
        df$cor[indx] <- "equal"
        df$max_rhat[indx] <- max(summary(model$model)$summary[, "Rhat"])
        df$min_neff[indx] <- min(summary(model$model)$summary[, "n_eff"])
        indx <- indx + 1
      }
    }

    if (length(which(variance %in% "unequal")) > 0) {
      for (i in seq(kmin, kmax)) {
        model <- fit_dfa(
          y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)),
          nu_fixed = 100, estimate_nu = FALSE, ...
        )
        df$num_trends[indx] <- i
        df$looic[indx] <- loo::loo(loo::extract_log_lik(model$model))$looic

        df$converge[indx] <- is_converged(model, convergence_threshold)
        # if model is best, keep it
        if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
          best_model <- model
          best_loo <- df$looic[indx]
        }
        df$error[indx] <- "normal"
        df$cor[indx] <- "independent"
        df$max_rhat[indx] <- max(summary(model$model)$summary[, "Rhat"])
        df$min_neff[indx] <- min(summary(model$model)$summary[, "n_eff"])
        indx <- indx + 1
      }
    }
  }

  df <- dplyr::arrange(df, .data$looic)

  # return best model = one that converges
  list(summary = df, best_model = best_model)
}
