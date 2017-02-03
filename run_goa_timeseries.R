library(rstan)
library(MARSS)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
source("find_dfa_trends.r")
# 
ak_full = read.csv("/users/eric.ward/downloads/fate data Litzow 2-1_working.csv")
ak_full = ak_full[which(!is.na(ak_full$year)),]
meta = read.csv("/users/eric.ward/downloads/fate data Litzow 2-1_meta.csv")
meta = meta[which(meta$System != ""),]

meta = meta[which(meta$System=="GOA"),]

ak_full = ak_full[,which(names(ak_full)%in%meta$Code)]

for(i in 1:ncol(ak_full)) ak_full[,i] = as.numeric(ak_full[,i])

df = find_dfa_trends(y = t(ak_full), iter=2000)

mod1 = fit_dfa(y = t(ak_full), num_trends = 2, iter=2000)
rt = rotate_trends(mod1)


par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(1950:2016, apply(rt$trends[,1,],2,mean), type="l", ylim=c(-4,4), lwd=3)
lines(1950:2016, apply(rt$trends[,1,],2,quantile,0.025,col="grey"))
lines(1950:2016, apply(rt$trends[,1,],2,quantile,0.975,col="grey"))
plot(1950:2016, apply(rt$trends[,2,],2,mean), type="l", ylim=c(-4,4), lwd=3)
lines(1950:2016, apply(rt$trends[,2,],2,quantile,0.025,col="grey"))
lines(1950:2016, apply(rt$trends[,2,],2,quantile,0.975,col="grey"))
plot(1950:2016, apply(rt$trends[,3,],2,mean), type="l", ylim=c(-4,4), lwd=3)
lines(1950:2016, apply(rt$trends[,3,],2,quantile,0.025,col="grey"))
lines(1950:2016, apply(rt$trends[,3,],2,quantile,0.975,col="grey"))

par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(1950:2015, diff(apply(rt$trends[,1,],2,mean)), type="l", lwd=3, ylab="anomaly")
plot(1950:2015, diff(apply(rt$trends[,2,],2,mean)), type="l", lwd=3, ylab="anomaly")
plot(1950:2015, diff(apply(rt$trends[,3,],2,mean)), type="l", lwd=3, ylab="anomaly")
xx = cbind(diff(apply(rt$trends[,1,],2,mean)), diff(apply(rt$trends[,2,],2,mean)), diff(apply(rt$trends[,3,],2,mean)))
plot(1950:2015, log(mvtnorm::dmvnorm(xx, mean = rep(0,3))), type="l", ylab="log(dmvnorm)")

par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(rollapply(zoo(apply(rt$trends[,1,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")
plot(rollapply(zoo(apply(rt$trends[,2,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")
plot(rollapply(zoo(apply(rt$trends[,3,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")

Z = rotate_trends(mod8)$Z_rot

Zmean = cbind(colnames(ak_full), round(apply(Z, c(2,3), mean), 2))

ggplot(df, aes(year, mean)) + 
  geom_ribbon(aes(ymin = low, ymax=high),alpha=0.4) + 
  facet_wrap(~trend) + geom_line()


# 2-stage model -- let's assume there's 3 trends for data
dfa1 = fit_dfa(y = t(ak_full[1:32,]), num_trends = 3, varIndx = seq(1,ncol(ak_full)))
zz_mean = apply(extract(dfa1)$Z, c(2,3), mean)
zz_sd = apply(extract(dfa1)$Z, c(2,3), sd)

z_mean = zz_mean[which(mat_indx != 0)]
z_sd = zz_sd[which(mat_indx != 0)]

data_list = list(
  N,
  P,
  K,
  nZ,
  y,
  row_indx,
  col_indx,
  nZero,
  row_indx_z,
  col_indx_z,
  nZero,
  row_indx_z,
  col_indx_z,
  row_indx_pos,
  col_indx_pos,
  n_pos,
  z_mean,
  z_sd
)

mod = stan(
  data = data_list,
  pars = c("x", "Z", "sigma", "log_lik"),
  file = "dfa_prior.stan",
  chains = chains,
  iter = iter,
  control = control
)

zz2_mean = apply(extract(mod)$Z, c(2,3), mean)
zz2_sd = apply(extract(mod)$Z, c(2,3), sd)

cor(zz_mean[,1], zz2_mean[,1])
cor(zz_sd[,1], zz2_sd[,1])

N = ncol(y)
P = nrow(y)
K = num_trends # number of dfa trends
nZ = P * K - sum(1:K) + K
mat_indx = matrix(0, P, 3)
start = 1
for (k in 1:3) {
  if (k == 1)
    mat_indx[, k] = (1:nZ)[start:(start + P - k)]
  if (k > 1)
    mat_indx[-c(0:(k - 1)), k] = (1:nZ)[start:(start + P - k)]
  start = start + (P - k + 1)
}

