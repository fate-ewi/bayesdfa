# Run example simulation

if (interactive()) options(mc.cores = parallel::detectCores())

set.seed(42)
num_trends <- 2
num_ts <- 5
num_years <- 25
loadings_matrix <- matrix(nrow = num_ts, ncol = num_trends,
  rnorm(num_ts * num_trends, 0, 0.6))
# loadings_matrix[loadings_matrix > 1] <- 1
# loadings_matrix[loadings_matrix > -1] <- -1

dat <- sim_dfa(
  num_trends = num_trends,
  num_years = num_years,
  num_ts = num_ts,
  loadings_matrix = loadings_matrix,
  sigma = rlnorm(1, meanlog = log(0.4), 0.4))

m <- fit_dfa(dat$y_sim, num_trends = num_trends, zscore = FALSE)

s <- reshape_samples(m$samples)
Zhat_m <- apply(s$Z, c(2, 3), median)
Zhat_l <- apply(s$Z, c(2, 3), quantile, probs = 0.025)
Zhat_u <- apply(s$Z, c(2, 3), quantile, probs = 0.975)
xhat_m <- apply(s$x, c(2, 3), median)
xhat_l <- apply(s$x, c(2, 3), quantile, probs = 0.025)
xhat_u <- apply(s$x, c(2, 3), quantile, probs = 0.975)

# Note I am assuming a single shared sigma here:
sigmahat_m <- apply(extract(m$model)$sigma, 2, median)
sigmahat_l <- apply(extract(m$model)$sigma, 2, quantile, probs = 0.025)
sigmahat_u <- apply(extract(m$model)$sigma, 2, quantile, probs = 0.975)

# ------------------------------
# Plot the trends:
df <- data.frame(
  x = c(t(xhat_m)),
  lo = c(t(xhat_l)),
  hi = c(t(xhat_u)),
  trend = paste0("Trend ", sort(rep(seq_len(num_trends), num_years))),
  time = rep(seq_len(num_years), num_trends),
  true = c(t(dat$x)))
p1 <- ggplot(df, aes_string(x = "time", y = "x")) +
  geom_ribbon(aes_string(ymin = "lo", ymax = "hi"), alpha = 0.4) +
  geom_line() + facet_wrap("trend") +
  xlab("Time") + ylab("") +
  geom_line(aes(y = true), col = "red")

# ------------------------------
# Plot the loadings:
df2 <- data.frame(
  z = c(Zhat_m),
  lo = c(Zhat_l),
  hi = c(Zhat_u),
  true = c(dat$Z),
  trend = rep(seq_len(num_trends), each = num_ts),
  ts = rep(seq_len(num_ts), num_trends))
p2 <- ggplot(df2, aes_string(x = "ts", y = "z")) +
  geom_pointrange(aes_string(ymin = "lo", ymax = "hi")) +
  facet_wrap("trend") +
  geom_point(aes(x = ts, y = true), col = "red", pch = 4, cex = 4) +
  xlab("Timeseries") +
  ylab("Loading")

gridExtra::grid.arrange(p1, p2)
