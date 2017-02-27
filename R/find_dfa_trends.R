#' Find the best number of trends according to LOO
#'
#' Fit a DFA with different number of trends and return the leave one out (LOO)
#' value as calculated by the loo package. Note that the models will be fit once
#' with shared variances and once with independent variances.
#'
#' @param y A matrix of data to fit. Columns represent time element.
#' @param kmin Minimum number of trends.
#' @param kmax Maximum number of trans.
#' @param iter Iterations when sampling from each Stan model.
#' @param compare_normal If TRUE, does model selection comparison of Normal vs
#'   Student-t errors
#' @param convergence_threshold The maximum allowed value of Rhat to determine
#'   convergence of parameters
#' @param variance Vector of variance arguments for searching over large groups of models. Can be either or both of ("equal","unequal")
#' @export
#'
#' @importFrom loo loo extract_log_lik
#' @importFrom stats quantile time varimax

find_dfa_trends = function(y = y, kmin = 1, kmax = 5, iter = 2000, compare_normal=TRUE, convergence_threshold=1.05, variance=c("equal","unequal")) {

  df = data.frame(
    model = seq(1, ifelse(compare_normal==FALSE, length(variance) * length(kmin:kmax), 2 * length(variance) * length(kmin:kmax))),
    num_trends = NA,
    looic = NA,
    cor = NA,
    error=NA,
    converge=FALSE,
    stringsAsFactors=FALSE
  )
  best_model = NULL
  best_loo = 1.0e50

  indx = 1
  
  if(length(which(variance%in%"equal")) > 0) {
  for (i in kmin:kmax) {
    model = fit_dfa(y = y, num_trends = i, iter = iter, nu_fixed = 7)
    df$num_trends[indx] = i
    df$looic[indx] = loo(extract_log_lik(model))$looic

    df$converge[indx] = is_converged(model, convergence_threshold)
    # if model is best, keep it
    if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
      best_model = model
      best_loo = df$looic[indx]
    }
    df$error[indx] = "student-t"
    df$cor[indx] = "equal"
    indx = indx + 1
  }
  }
  
  if(length(which(variance%in%"unequal")) > 0) {
  for (i in kmin:kmax) {
    model = fit_dfa(y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)), nu_fixed=7
    )
    df$num_trends[indx] = i
    df$looic[indx] = loo::loo(loo::extract_log_lik(model))$looic

    df$converge[indx] = is_converged(model, convergence_threshold)
    # if model is best, keep it
    if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
      best_model = model
      best_loo = df$looic[indx]
    }
    df$error[indx] = "student-t"
    df$cor[indx] = "independent"
    indx = indx + 1
  }
  }
  
  if(length(which(variance%in%"equal")) > 0) {
  if(compare_normal==TRUE) {
    for (i in kmin:kmax) {
      model = fit_dfa(y = y, num_trends = i, iter = iter, nu_fixed = 100)
      df$num_trends[indx] = i
      df$looic[indx] = loo(extract_log_lik(model))$looic

      df$converge[indx] = is_converged(model, convergence_threshold)
      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
        best_model = model
        best_loo = df$looic[indx]
      }
      df$error[indx] = "normal"
      df$cor[indx] = "equal"
      indx = indx + 1
    }
    }
    
    if(length(which(variance%in%"unequal")) > 0) {
    for (i in kmin:kmax) {
      model = fit_dfa(y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)), nu_fixed=100
      )
      df$num_trends[indx] = i
      df$looic[indx] = loo::loo(loo::extract_log_lik(model))$looic

      df$converge[indx] = is_converged(model, convergence_threshold)
      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == TRUE) {
        best_model = model
        best_loo = df$looic[indx]
      }
      df$error[indx] = "normal"
      df$cor[indx] = "independent"
      indx = indx + 1
      }
     }
   }  

  df <- dplyr::arrange_(df, ~ looic)

  # return best model = one that converges
  list(summary = df, best_model = best_model)
}
