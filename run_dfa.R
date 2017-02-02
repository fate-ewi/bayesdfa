library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
data(harborSealWA)
y = t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

fit1 = fit_dfa(y = y, num_trends = 1)
fit2 = fit_dfa(y = y, num_trends = 2)
fit3 = fit_dfa(y = y, num_trends = 3)

# Illustrate how to get trends 
rotated1 = rotate_trends(fit2)
matplot(rotated1$trends)

