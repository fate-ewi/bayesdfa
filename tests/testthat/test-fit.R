library(MARSS)
if (interactive()) {
  options(mc.cores = parallel::detectCores())
}

data(harborSealWA)
y <- t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

fit1 <- fit_dfa(y = y, num_trends = 1)

testthat::expect_output(print(fit1), "Inference for Stan model")
