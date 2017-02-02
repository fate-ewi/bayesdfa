find_dfa_trends = function(y = y, k = 1:5, iter=2000) {

  df = data.frame("model"=2*length(k), "num_trends"=NA, "looic"=NA, "cor"=NA)
  indx = 1
  for(i in min(k):max(k)) {
  model = fit_dfa(y = y, num_trends = k, iter=iter)
  df$num_trends[i] = i
  df$looic[i] = loo(extract_log_lik(model))$looic
  df$cor[i] = "equal"
  indx = indx + 1
  }

  for(i in min(k):max(k)) {
    model = fit_dfa(y = y, num_trends = k, iter=iter, varIndx = seq(1,nrow(y)))
    df$num_trends[i] = i
    df$looic[i] = loo(extract_log_lik(model))$looic
    df$cor[i] = "independent"
    indx = indx + 1
  }
  return(df)
}