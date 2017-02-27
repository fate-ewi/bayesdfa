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
#' @param compare_normal If TRUE, does model selection comparison of Normal vs Student-t errors
#' @export
#'
#' @importFrom loo loo extract_log_lik
#' @importFrom stats quantile time varimax

find_dfa_trends = function(y = y, kmin = 1, kmax = 5, iter = 2000, compare_normal=TRUE) {

  df = data.frame(
    model = seq(1, ifelse(compare_normal==FALSE, 2 * length(kmin:kmax), 4 * length(kmin:kmax))),
    num_trends = NA,
    looic = NA,
    cor = NA,
    error=NA,
    converge="N",
    stringsAsFactors=FALSE
  )
  best_model = NULL
  best_loo = 1.0e50

  indx = 1
  for (i in kmin:kmax) {
    model = fit_dfa(y = y, num_trends = i, iter = iter, nu = 7)
    df$num_trends[indx] = i
    df$looic[indx] = loo(extract_log_lik(model))$looic

    Rhats = summary(model)$summary[,"Rhat"]
    if(max(Rhats[grep("Z",names(Rhats))], na.rm=T) < 1.05) df$converge[indx] = "Y"
    # if model is best, keep it
    if (df$looic[indx] < best_loo & df$converge[indx] == "Y") {
      best_model = model
      best_loo = df$looic[indx]
    }
    df$error[indx] = "student-t"
    df$cor[indx] = "equal"
    indx = indx + 1
  }

  for (i in kmin:kmax) {
    model = fit_dfa(y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)), nu=7
    )
    df$num_trends[indx] = i
    df$looic[indx] = loo::loo(loo::extract_log_lik(model))$looic

    Rhats = summary(model)$summary[,"Rhat"]
    if(max(Rhats[grep("Z",names(Rhats))], na.rm=T) < 1.05) df$converge[indx] = "Y"
    # if model is best, keep it
    if (df$looic[indx] < best_loo & df$converge[indx] == "Y") {
      best_model = model
      best_loo = df$looic[indx]
    }
    df$error[indx] = "student-t"
    df$cor[indx] = "independent"
    indx = indx + 1
  }

  if(compare_normal==TRUE) {
    for (i in kmin:kmax) {
      model = fit_dfa(y = y, num_trends = i, iter = iter, nu = 100)
      df$num_trends[indx] = i
      df$looic[indx] = loo(extract_log_lik(model))$looic
      
      Rhats = summary(model)$summary[,"Rhat"]
      if(max(Rhats[grep("Z",names(Rhats))], na.rm=T) < 1.05) df$converge[indx] = "Y"
      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == "Y") {
        best_model = model
        best_loo = df$looic[indx]
      }
      df$error[indx] = "normal"
      df$cor[indx] = "equal"
      indx = indx + 1
    }
    
    for (i in kmin:kmax) {
      model = fit_dfa(y = y, num_trends = i, iter = iter, varIndx = seq(1, nrow(y)), nu=100
      )
      df$num_trends[indx] = i
      df$looic[indx] = loo::loo(loo::extract_log_lik(model))$looic
      
      Rhats = summary(model)$summary[,"Rhat"]
      if(max(Rhats[grep("Z",names(Rhats))], na.rm=T) < 1.05) df$converge[indx] = "Y"
      # if model is best, keep it
      if (df$looic[indx] < best_loo & df$converge[indx] == "Y") {
        best_model = model
        best_loo = df$looic[indx]
      }
      df$error[indx] = "normal"
      df$cor[indx] = "independent"
      indx = indx + 1
    }
  }
  
  df <- dplyr::arrange_(df, ~ looic)

  # return best model = one that converges
  list(summary = df, best_model = best_model)
}
