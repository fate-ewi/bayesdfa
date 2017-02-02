library(rstan)
library(MARSS)
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

df = find_dfa_trends(y = t(ak_full))
# 5 models with equal errors
mod1 = fit_dfa(y = t(ak_full), num_trends = 1, iter=2000)
mod2 = fit_dfa(y = t(ak_full), num_trends = 2, iter=2000)
mod3 = fit_dfa(y = t(ak_full), num_trends = 3, iter=2000)
mod4 = fit_dfa(y = t(ak_full), num_trends = 4, iter=2000)
mod5 = fit_dfa(y = t(ak_full), num_trends = 5, iter=2000)

# same 5 models with independent errors
mod6 = fit_dfa(y = t(ak_full), num_trends = 1, iter=2000, varIndx = seq(1,ncol(ak_full)))
mod7 = fit_dfa(y = t(ak_full), num_trends = 2, iter=2000, varIndx = seq(1,ncol(ak_full)))
mod8 = fit_dfa(y = t(ak_full), num_trends = 3, iter=2000, varIndx = seq(1,ncol(ak_full)))
mod9 = fit_dfa(y = t(ak_full), num_trends = 4, iter=2000, varIndx = seq(1,ncol(ak_full)))
mod10 = fit_dfa(y = t(ak_full), num_trends = 5, iter=2000, varIndx = seq(1,ncol(ak_full)))

# best model = independent errors, 3 trends
rt = rotate_trends(mod8)  

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



