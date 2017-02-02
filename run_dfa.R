library(rstan)
library(MARSS)

source("fit_dfa.r")

fit1 = fit_dfa(y =y, num_trends = 1)
fit2 = fit_dfa(y =y, num_trends = 2)
fit3 = fit_dfa(y =y, num_trends = 3)


