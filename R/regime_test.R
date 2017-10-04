library(rstan)
data(Nile)
y = log(Nile)

pars = c("u", "pred", "sigma", "log_lik")

# fit null model with 1 regime
data_list = list(N = length(y), y = y)
fit_1 = stan("exec/regime_1.stan", data = data_list, iter = 2000,
  chains = 3, verbose = TRUE)
loo::loo(loo::extract_log_lik(fit_1))$looic

# fit model with 2 regimes
data_list = list(N = length(y), y = y)
fit_2 = stan("exec/regime_2.stan", data = data_list, iter = 2000,
  chains = 3, verbose = TRUE)
loo::loo(loo::extract_log_lik(fit_2))$looic

# fit the model with 3 regimes
regimes=3
data_list = list(N = length(y), n_regime = regimes, y = y, ones = rep(1, regimes-1))
fit_3 = stan("exec/regime_3plus.stan", data = data_list, iter = 2000,
  chains = 3, verbose = TRUE)
loo::loo(loo::extract_log_lik(fit_3))$looic

regimes=4
data_list = list(N = length(y), n_regime = regimes, y = y, ones = rep(1, regimes-1))
fit_4 = stan("exec/regime_3plus.stan", data = data_list, iter = 2000,
  chains = 3, verbose = TRUE)
loo::loo(loo::extract_log_lik(fit_4))$looic
