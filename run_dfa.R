library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
data(harborSealWA)
y = t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

fit1 = fit_dfa(y = y, num_trends = 1)
rstan::traceplot(fit1, pars = "Z")
fit2 = fit_dfa(y = y, num_trends = 2)
rstan::traceplot(fit2, pars = "Z")
fit3 = fit_dfa(y = y, num_trends = 3)
rstan::traceplot(fit3, pars = "Z")

# Illustrate how to get trends 
rotated1 = rotate_trends(fit2)
matplot(rotated1$trends)

