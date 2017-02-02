find_dfa_trends = function(y = y, kmin=1, kmax=5, iter=2000) {

  df = data.frame("model"=seq(1,2*length(kmin:kmax)), "num_trends"=NA, "looic"=NA, "cor"=NA)
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