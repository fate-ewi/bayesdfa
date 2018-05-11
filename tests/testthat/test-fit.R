context("Fitting")

if (interactive()) options(mc.cores = parallel::detectCores())

y <- t(scale(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")]))

set.seed(1)
fit1 <- fit_dfa(y = y, num_trends = 1, iter = 800, chains = 1)

test_that("MARSS and bayesdfa match", {
  ml_fit <- MARSS::MARSS(y, form = "dfa", model = list(m = 1))
  ml_means <- c(ml_fit$states)
  bayes_means <- apply(extract(fit1$model, "x")$x[, 1, ], 2, mean)
  expect_equal(cor(abs(bayes_means), abs(ml_means)), 1, tolerance = 0.01)
})

test_that("print method works", {
  expect_output(print(fit1), "n_eff")
})

test_that("NA indexing works", {
  yy <- matrix(nrow = 3, ncol = 3, data = 1)
  yy[1, 1] <- NA
  yy[2, 3] <- NA
  m <- fit_dfa(yy, num_trends = 1, sample = FALSE, zscore = FALSE)
  expect_equal(m$n_na, 2L)
  expect_equal(m$row_indx_na, c(1L, 2L))
  expect_equal(m$col_indx_na, c(1L, 3L))
})

test_that("if time series all same value then zscore stops with error", {
  yy <- matrix(nrow = 3, ncol = 3, data = 1)
  expect_error(fit_dfa(y = yy, num_trends = 1, sample = FALSE, zscore = TRUE))
})

test_that("find_dfa_trends works", {
  skip_on_cran()

  set.seed(42)
  s <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 3)
  x <- find_dfa_trends(
    y = s$y_sim, iter = 1000,
    kmin = 1, kmax = 2, chains = 1, compare_normal = FALSE,
    variance = "equal", convergence_threshold = 1.1,
    control = list(adapt_delta = 0.95, max_treedepth = 20)
  )

  expect_equal(x$summary$model, c(2L, 1L))
  expect_lt(x$summary$looic[[1]], x$summary$looic[[2]])
})

# test_that("we can find which chains to flip", {
#   skip_on_cran()
#   skip_on_travis()
#   skip_on_appveyor()
#
# set.seed(1)
# m <- fit_dfa(y = y, num_trends = 3, iter = 1000, chains = 4)
# expect_equal(find_inverted_chains(m, trend = 1, plot = TRUE), 4)
#
# # a single chain:
# m <- fit_dfa(y = y, num_trends = 3, iter = 1000, chains = 1)
# expect_equal(find_inverted_chains(m, trend = 1, plot = TRUE), 1)
# })
