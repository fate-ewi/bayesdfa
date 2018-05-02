if (interactive()) options(mc.cores = parallel::detectCores())

# set.seed(1)
# num_trends <- 2
# num_ts <- 3
# num_years <- 30
#
#
# dat <- sim_dfa(
#   num_trends = num_trends,
#   num_years = num_years,
#   num_ts = num_ts,
#   loadings_matrix = loadings_matrix,
#   sigma = 0.2, nu_fixed = 200)
#
# m2 <- fit_dfa(dat$y_sim, num_trends = num_trends, zscore = TRUE,
#   iter = 1000, chains = 4, seed = 1)
#
# x <- rotate_trends(m2)
# plot_trends(x)
# plot_loadings(x)
