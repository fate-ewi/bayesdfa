if (interactive()) options(mc.cores = parallel::detectCores())

y <- t(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

set.seed(1)

test_that("Basic model fits", {
  fit1 <- fit_dfa(y = y, num_trends = 1, iter = 1000, chains = 1)
  expect_output(print(fit1), "n_eff")
})

# test_that("find_dfa_trends works", {
#   skip_on_cran()
#   skip_on_travis()
#   skip_on_appveyor()
#   dfa_summary <- find_dfa_trends(y = y, iter = 1000, kmin = 1, kmax = 1)
# })

test_that("we can find which chains to flip", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()

  set.seed(1)
  # m <- fit_dfa(y = y, num_trends = 3, iter = 1000, chains = 4)
  # expect_equal(find_inverted_chains(m, trend = 1, plot = TRUE), 4)
  #
  # # a single chain:
  # m <- fit_dfa(y = y, num_trends = 3, iter = 1000, chains = 1)
  # expect_equal(find_inverted_chains(m, trend = 1, plot = TRUE), 1)
})
