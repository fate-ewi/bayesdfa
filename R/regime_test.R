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

e <- extract(m)
ll <- plyr::laply(1:500, function(i) {
  dnorm(y, mean = c(e$e[seq(1, e$s[i])-1], e$l[seq(e$s[i], length(y))]),
    sd = e$sigma[i], log = TRUE)
})
loo::loo(ll)

##########

y <- c(rnorm(10, 1.5, 0.2), rnorm(60, 1, 0.2), rnorm(10, 1.5, 0.2))
plot(y)
initf <- function() list(e = 2, l = 0, sigma = 0.2)
m <- stan("R/cp_norm.stan",
  data = list(T = length(y), D = y),
  cores = 1, chains = 1, iter = 1000,
  init = initf)

hist(extract(m)$s)
x <- extract(m)$expectations
plot(apply(x, 2, mean), type = "l")
hist(extract(m)$l)
hist(extract(m)$e)
hist(extract(m)$sigma)
e <- extract(m)
plot(e$l, e$s)
plot(e$e, e$s)
print(m, pars = c("e", "l", "sigma"))
######

set.seed(1234)
y <- c(rnorm(15, 1.2, 0.3), rnorm(15, -0.5, 0.3), rnorm(15, 2, 0.3))
plot(y)
m <- stan("R/cp_norm_multiple.stan",
  data = list(T = length(y), D = y),
  cores = 3, chains = 3, iter = 300,
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

plot(y)
plot(apply(x, 2, mean), type = "l")
abline(v = 16)
abline(v = 31)

# sigma:
y <- c(rnorm(25, 1.2, 0.3), rnorm(25, -0.5, 0.3))
plot(y)
m <- stan("R/cp_norm_sigma.stan",
  data = list(T = length(y), D = y, sigma = 0.3),
  cores = 8, chains = 8, iter = 500)
hist(extract(m)$s)
x <- extract(m)$expectations
plot(apply(x, 2, mean), type = "l")
hist(extract(m)$l)
hist(extract(m)$e)
hist(extract(m)$sigma)

# t with sigma
x <- rnorm(100, qnorm(0.025), 0.2)
hist(x, xlim = c(-3, 3))
p <- pnorm(x, mean = 0, sd = 1)
hist(p)
round(median(p), 3)

