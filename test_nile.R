library(bayesdfa)
data(Nile)
y = log(Nile)

# extreme demo
y[c(1:10, 33:40, 62:75)] = y[c(1:10, 33:40, 62:75)] + 3

fit = find_regimes(y = y, n_regimes = 1)
loo::loo(loo::extract_log_lik(fit$model))$looic

fit = find_regimes(y = y, n_regimes = 2)
loo::loo(loo::extract_log_lik(fit$model))$looic

fit = find_regimes(y = y, n_regimes = 3)
loo::loo(loo::extract_log_lik(fit$model))$looic

fit = find_regimes(y = y, n_regimes = 4)
loo::loo(loo::extract_log_lik(fit$model))$looic

#fit = find_regimes(y = y, n_regimes = 5)
#loo::loo(loo::extract_dlog_lik(fit$model))$looic

#fit = find_regimes(y = y, n_regimes = 6)
#loo::loo(loo::extract_log_lik(fit$model))$looic
