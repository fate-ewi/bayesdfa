library(MARSS)
data(harborSealWA)
y <- t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

set.seed(1)

test_that("Basic model fits", {
  fit1 <- fit_dfa(y = y, num_trends = 1, iter = 1000)
  expect_output(print(fit1), "Inference for Stan model")
})

test_that("find_dfa_trends works", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  dfa_summary <- find_dfa_trends(y = y, iter = 1000, kmin = 1, kmax = 1)
})
