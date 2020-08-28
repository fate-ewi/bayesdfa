#' Find the best number of trends according to LOOIC
#'
#' Fit a DFA with different number of trends and return the leave one out (LOO)
#' value as calculated by the [loo][loo::loo-package()] package.
#'
#' @param y A matrix of data to fit. Columns represent time element.
#' @param kmin Minimum number of trends, defaults to 1.
#' @param kmax Maximum number of trends, defaults to 5.
#' @param iter Iterations when sampling from each Stan model, defaults to 2000.
#' @param thin Thinning rate when sampling from each Stan model, defaults to 1.
#' @param compare_normal If `TRUE`, does model selection comparison of Normal vs.
#'   Student-t errors
#' @param convergence_threshold The maximum allowed value of Rhat to determine
#'   convergence of parameters
#' @param variance Vector of variance arguments for searching over large groups
#'   of models. Can be either or both of ("equal","unequal")
#' @param ... Other arguments to pass to `fit_dfa()`
#' @export
#' @examples
#' \donttest{
#' set.seed(42)
#' s <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 3)
#' # only 1 chain and 180 iterations used so example runs quickly:
#' m <- find_dfa_trends(
#'   y = s$y_sim, iter = 50,
#'   kmin = 1, kmax = 2, chains = 1, compare_normal = FALSE,
#'   variance = "equal", convergence_threshold = 1.1,
#'   control = list(adapt_delta = 0.95, max_treedepth = 20))
#' m$summary
#' m$best_model
#' }
#' @importFrom loo loo extract_log_lik
#' @importFrom stats quantile time varimax
#' @importFrom rlang .data

find_dfa_trends <- function(y = y,
  kmin = 1,
  kmax = 5,
  iter = 2000,
  thin = 1,
  compare_normal = FALSE,
  convergence_threshold = 1.05,
  variance = c("equal", "unequal"),
  ...) {

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
        y = y, num_trends = i, iter = iter, thin = thin,
        estimate_nu = TRUE, ...
      )

      df$converge[indx] <- is_converged(model, convergence_threshold)
      df$num_trends[indx] <- i

      # relative effective sample size
      log_lik <- loo::extract_log_lik(model$model, merge_chains = FALSE)
      n_chains <- dim(rstan::extract(model$model, "log_lik", permuted=FALSE))[2]
      rel_eff <- loo::relative_eff(exp(log_lik))
      # calculate looic
      df$looic[indx] <- loo::loo(log_lik, r_eff = rel_eff)$estimates["looic",1]

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
        y = y, num_trends = i, iter = iter, thin = thin, varIndx = seq(1, nrow(y)),
        estimate_nu = TRUE, ...
      )
      df$num_trends[indx] <- i

      log_lik <- loo::extract_log_lik(model$model, merge_chains = FALSE)
      n_chains <- dim(rstan::extract(model$model, "log_lik", permuted=FALSE))[2]
      rel_eff <- loo::relative_eff(exp(log_lik))
      # calculate looic
      df$looic[indx] <- loo::loo(log_lik, r_eff = rel_eff)$estimates["looic",1]

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
          y = y, num_trends = i, iter = iter, thin = thin, nu_fixed = 100,
          estimate_nu = FALSE, ...
        )
        df$num_trends[indx] <- i

        log_lik <- loo::extract_log_lik(model$model, merge_chains = FALSE)
        n_chains <- dim(rstan::extract(model$model, "log_lik", permuted=FALSE))[2]
        rel_eff <- loo::relative_eff(exp(log_lik))
        # calculate looic
        df$looic[indx] <- loo::loo(log_lik, r_eff = rel_eff)$estimates["looic",1]

        df$converge[indx] <- is_converged(model, convergence_threshold)
        # if model is best, keep it
        if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
          best_model <- model
          best_loo <- df$looic[indx]
        }
        df$error[indx] <- "normal"
        df$cor[indx] <- "equal"
        #df$max_rhat[indx] <- max(as.data.frame(summary(model$model)$summary)[,"Rhat"])
        #df$min_neff[indx] <- min(as.data.frame(summary(model$model)$summary)[,"n_eff"])
        indx <- indx + 1
      }
    }

    if (length(which(variance %in% "unequal")) > 0) {
      for (i in seq(kmin, kmax)) {
        model <- fit_dfa(
          y = y, num_trends = i, iter = iter, thin = thin, varIndx = seq(1, nrow(y)),
          nu_fixed = 100, estimate_nu = FALSE, ...
        )
        df$num_trends[indx] <- i

        log_lik <- loo::extract_log_lik(model$model, merge_chains = FALSE)
        n_chains <- dim(rstan::extract(model$model, "log_lik", permuted=FALSE))[2]
        rel_eff <- loo::relative_eff(exp(log_lik))
        # calculate looic
        df$looic[indx] <- loo::loo(log_lik, r_eff = rel_eff)$estimates["looic",1]

        df$converge[indx] <- is_converged(model, convergence_threshold)
        # if model is best, keep it
        if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
          best_model <- model
          best_loo <- df$looic[indx]
        }
        df$error[indx] <- "normal"
        df$cor[indx] <- "independent"
        #df$max_rhat[indx] <- max(as.data.frame(summary(model$model)$summary)[,"Rhat"])
        #df$min_neff[indx] <- min(as.data.frame(summary(model$model)$summary)[,"n_eff"])
        indx <- indx + 1
      }
    }
  }

  df <- dplyr::arrange(df, .data$looic)

  # return best model = one that converges
  list(summary = df, best_model = best_model)
}
