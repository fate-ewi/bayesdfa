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

ggplot(df, aes(year, mean)) + 
  geom_ribbon(aes(ymin = low, ymax=high),alpha=0.4) + 
  facet_wrap(~trend) + geom_line()

meta = meta[which(meta$System=="GOA"),]
ak_full = ak_full[,which(names(ak_full)%in%meta$Code)]
ak_full = ak_full[,-which(names(ak_full)== "ICY.ZOOP.TOTDEN")]
for(i in 1:ncol(ak_full)) ak_full[,i] = as.numeric(ak_full[,i])

# Run 3-trend model dropping out batches of data from same survey
mod1 = fit_dfa(t(ak_full[,12:36]), num_trends = 3, varIndx = seq(1,25)) # no salmon
mod2 = fit_dfa(t(ak_full[,c(1:11,18:36)]), num_trends = 3, varIndx = seq(1,30)) # no salmon
mod3 = fit_dfa(t(ak_full[,c(1:17,23:36)]), num_trends = 3, varIndx = seq(1,31)) # no salmon
mod4 = fit_dfa(t(ak_full[,c(1:31)]), num_trends = 3, varIndx = seq(1,31)) # no salmon

rt1 = rotate_trends(mod1)
rt2 = rotate_trends(mod2)
rt3 = rotate_trends(mod3)
rt4 = rotate_trends(mod4)

pdf("4 trends.pdf")
par(mfrow=c(2,2),mai=c(0.5,0.5,0.1,0.1))
matplot(t(rt1$trends_mean), main="No salmon")
matplot(t(rt2$trends_mean), main = "No Pav")
matplot(t(rt3$trends_mean), main = "No Chi")
matplot(t(rt4$trends_mean), main = "No Icy")
dev.off()