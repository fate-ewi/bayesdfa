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

set.seed(1234)
y <- c(rnorm(15, 1.2, 0.3), rnorm(15, 0.5, 0.3), rnorm(15, -0.5, 0.3))
plot(y)
m <- stan("R/cp_norm_multiple.stan",
  data = list(T = length(y), D = y),
  cores = 3, chains = 3, iter = 200,
  pars = "lp", include = FALSE)
# m
print(m, pars = c("e", "m", "l", "sigma"))
# hist(extract(m)$s)
x <- extract(m)$expectations
ex <- matrix(data = 0, nrow = dim(x)[2], ncol = dim(x)[3])
for (i in 1:dim(x)[2]) {
  for (j in 1:dim(x)[3]) {
    ex[i, j] <- mean(x[,i,j])
  }}

image(x = seq_along(y), y = seq_along(y), z = ex)
abline(v = 16)
abline(v = 31)

plot(apply(ex, 1, mean), type = "l")
