# library(rstan)
# x <- c(rep(0, 10), rep(1, 10))
# y <- rnorm(n = length(x), mean = 0.2 + x * 0.4, 0.5)
# plot(x, y)
# m <- rstan::sampling(stanmodels$ttest, data = list(x = x, y = y, N = length(x)),
#   iter = 400, chains = 1)
# m

#
# set.seed(11)
# y <- c(rep(1, 20), rep(4, 30), rep(1, 20)) + rnorm(70, 0, 0.5)
# plot(1:length(y), y)
# i <- window_minimum
# regime_years <- 0
#
# while (i < length(y)) {
#   # length_current_regime <- 0
#   possible_new_regime <- FALSE
#   new_regime <- FALSE
#   window_minimum <- 10
#   current_window <- 0
#   regime_tracker <- logical(length = length(y))
#
#   while (!new_regime) {
#     i <- i + 1
#     current_window <- current_window + 1
#     if (current_window >= window_minimum) {
#       x <- c(rep(0, 10), rep(1, current_window))
#       m <- lm(y[(regime_years[length(regime_years)] + 1):i] ~ x)
#       p <- summary(m)$coef[2, 4]
#       if (p < 0.05 & !possible_new_regime) {
#         possible_new_regime <- TRUE
#         regime_tracker[i] <- TRUE
#       }
#     }
#     if (possible_new_regime) {
#       new_regime <- TRUE
#       regime_years <- c(regime_years, i)
#     }
#   }
# }
# plot(1:length(y), y)
# abline(v = seq_along(regime_tracker)[regime_tracker])
#
