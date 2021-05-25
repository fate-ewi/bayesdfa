context("Fitting")

if (interactive()) options(mc.cores = parallel::detectCores())

set.seed(1)

# marss_installed <- "MARSS" %in% rownames(installed.packages())
# if (marss_installed) {
#   y <- t(scale(MARSS::harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")]))
#   fit1 <- fit_dfa(y = y, num_trends = 1, iter = 600, chains = 1)
#   test_that("MARSS and bayesdfa match", {
#     ml_fit <- MARSS::MARSS(y, form = "dfa", model = list(m = 1))
#     ml_means <- c(ml_fit$states)
#     bayes_means <- apply(extract(fit1$model, "x")$x[, 1, ], 2, mean)
#     expect_equal(cor(abs(bayes_means), abs(ml_means)), 1, tolerance = 0.01)
#   })
#
#   test_that("print method works", {
#     expect_output(print(fit1), "n_eff")
#   })
# }

test_that("est_correlation = TRUE works", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 3)
  m <- fit_dfa(
    y = s$y_sim, iter = 50, chains = 1,
    num_trends = 2, est_correlation = TRUE
  )
  expect_equal(class(m$model)[[1]], "stanfit")
})

test_that("NA indexing works", {
  yy <- matrix(nrow = 3, ncol = 3, data = 1)
  yy[1, 1] <- NA
  yy[2, 3] <- NA
  m <- fit_dfa(yy, num_trends = 1, sample = FALSE, scale = "center")
  expect_equal(m$sampling_args$data$n_na, 2L)
  expect_equal(m$sampling_args$data$row_indx_na, c(1L, 2L))
  expect_equal(m$sampling_args$data$col_indx_na, c(1L, 3L))
})

test_that("if time series all same value then zscore stops with error", {
  yy <- matrix(nrow = 3, ncol = 3, data = 1)
  expect_error(fit_dfa(y = yy, num_trends = 1, sample = FALSE, scale = "zscore"))
})

test_that("find_dfa_trends works", {
  skip_on_cran()

  set.seed(42)
  s <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 3)
  expect_warning({
    x <- find_dfa_trends(
      y = s$y_sim, iter = 200,
      kmin = 1, kmax = 2, chains = 1, compare_normal = FALSE,
      variance = "equal", convergence_threshold = 1.1,
      control = list(adapt_delta = 0.95, max_treedepth = 20)
    )
  })

  expect_equal(x$summary$model, c(2L, 1L))
  expect_lt(x$summary$looic[[1]], x$summary$looic[[2]])
})

test_that("long format data works", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
  m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 1, seed = 42)
  wide_means <- apply(extract(m$model, "x")$x[, 1, ], 2, mean)
  # fit long format data
  long <- data.frame(
    "obs" = c(s$y_sim[1, ], s$y_sim[2, ], s$y_sim[3, ]),
    "ts" = sort(rep(1:3, 20)), "time" = rep(1:20, 3)
  )
  m2 <- fit_dfa(y = long, data_shape = "long", iter = 100, chains = 1, num_trends = 1, seed = 42)
  long_means <- apply(extract(m2$model, "x")$x[, 1, ], 2, mean)
  expect_equal(cor(wide_means, long_means), 1, tolerance = 0.01)
})

test_that("compositional model works", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
  m <- fit_dfa(
    y = s$y_sim, iter = 50, chains = 1, num_trends = 2, seed = 42,
    z_model = "proportion"
  )

  expect_equal(class(m$model)[[1]], "stanfit")
})

test_that("compositional model works_2", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 2, num_years = 20, num_ts = 3)
  m <- fit_dfa(
    y = s$y_sim, iter = 50, chains = 1, num_trends = 2, seed = 42,
    z_model = "proportion"
  )

  expect_equal(class(m$model)[[1]], "stanfit")
})

test_that("estimate_sigma_process_1", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
  m <- fit_dfa(
    y = s$y_sim, iter = 50, chains = 1, num_trends = 2, seed = 42,
    estimate_process_sigma = TRUE, equal_process_sigma = TRUE
  )

  expect_equal(class(m$model)[[1]], "stanfit")
})

test_that("estimate_sigma_process_k", {
  skip_on_cran()
  set.seed(42)
  s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
  m <- fit_dfa(
    y = s$y_sim, iter = 50, chains = 1, num_trends = 2, seed = 42,
    estimate_process_sigma = TRUE, equal_process_sigma = FALSE
  )

  expect_equal(class(m$model)[[1]], "stanfit")
})
#
# test_that("estimate_spline_model", {
#   skip_on_cran()
#   set.seed(42)
#   s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#   m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 2, seed = 42,
#     estimate_process_sigma = TRUE, n_knots = 10, trend_model = "spline")
#
#   expect_equal(class(m$model)[[1]], "stanfit")
# })
#
# test_that("estimate_gp_model", {
#   skip_on_cran()
#   set.seed(42)
#   s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#   m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 2, seed = 42,
#     estimate_process_sigma = TRUE, n_knots = 5, trend_model = "gp")
#   expect_equal(class(m$model)[[1]], "stanfit")
# })
#
# test_that("estimate_rw_model_pars", {
#   skip_on_cran()
#   set.seed(42)
#   s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#   m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 1, seed = 42,
#     estimate_process_sigma = TRUE)
#   x_mean = apply(rstan::extract(m$model,"x")$x[,1,], 2, mean)
#   true_x = c(0.17, -0.40,  0.42, -1.45, -0.16, -1.85, -1.35, -0.62,
#     -2.15, -4.69, -5.22, -5.81, -2.81, -0.04,  3.35,  5.22,  7.55,  5.21,  2.51,  1.53)
#   expect_lt(sum(abs(true_x - x_mean)), 0.06)
# })
#
# test_that("estimate_gp_model_pars", {
#   skip_on_cran()
#   set.seed(42)
#   s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#   m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 2, seed = 42,
#     estimate_process_sigma = TRUE, n_knots = 5, trend_model = "gp")
#   x_mean = apply(rstan::extract(m$model,"x")$x[,1,], 2, mean)
#   true_x = c(-0.06, -0.08, -0.08, -0.10, -0.15, -0.29, -0.50, -0.84, -1.36, -2.16, -2.04,
#     -0.98, -0.13,  0.68,  1.63,  1.47,  1.11,  0.95,  0.93,  1.09)
#   expect_lt(sum(abs(true_x - x_mean)), 0.06)
# })
#
# test_that("estimate_spline_model_pars", {
#   skip_on_cran()
#   set.seed(42)
#   s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#   m <- fit_dfa(y = s$y_sim, iter = 100, chains = 1, num_trends = 1, seed = 42,
#     estimate_process_sigma = TRUE, n_knots = 10, trend_model = "spline")
#   x_mean = apply(rstan::extract(m$model,"x")$x[,1,], 2, mean)
#   true_x = c(0.34,  0.05, -0.17, -0.31, -0.37, -0.39, -0.42,
#     -0.58, -1.07, -2.00, -2.86, -2.89, -1.71,  0.12,  1.87,  3.01,  3.27,  2.46,  1.14,  0.30)
#   expect_lt(sum(abs(true_x - x_mean)), 0.06)
# })

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
