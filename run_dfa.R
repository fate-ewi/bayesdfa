library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")

data(harborSealWA)
y = t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

fit1 = fit_dfa(y = y, num_trends = 1)
fit2 = fit_dfa(y = y, num_trends = 2)
fit3 = fit_dfa(y = y, num_trends = 3)
