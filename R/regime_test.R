library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

y <- c(rpois(50, 2.5), rpois(50, 0.5))
plot(y)
m <- stan("R/cp.stan",
  data = list(r_e = 2.5, r_l = 0.5, T = length(y), D = y),
  cores = 2, chains = 2, iter = 500)
hist(extract(m)$s)
x <- extract(m)$expectations
plot(apply(x, 2, mean), type = "l")


y <- c(rnorm(50, 1.5, 0.3), rnorm(50, 0.5, 0.3))
plot(y)
m <- stan("R/cp_norm.stan",
  data = list(T = length(y), D = y),
  cores = 2, chains = 2, iter = 500)
hist(extract(m)$s)
x <- extract(m)$expectations
plot(apply(x, 2, mean), type = "l")
hist(extract(m)$l)
hist(extract(m)$e)
hist(extract(m)$sigma)

set.seed(1)
y <- c(rnorm(15, 1, 0.3), rnorm(15, 0.5, 0.3), rnorm(15, -.5, 0.3))
plot(y)
m <- stan("R/cp_norm_multiple.stan",
  data = list(T = length(y), D = y),
  cores = 4, chains = 4, iter = 250,
  pars = "lp", include = FALSE)
print(m)
# hist(extract(m)$s)
# x <- extract(m)$expectations
# plot(apply(x, 2, mean), type = "l")
