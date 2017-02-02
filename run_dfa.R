library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
data(harborSealWA)
y = t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

fit1 = fit_dfa(y = y, num_trends = 1)
rstan::traceplot(fit1, pars = "Z")
fit2 = fit_dfa(y = y, num_trends = 2)
rstan::traceplot(fit2, pars = "Z")
fit3 = fit_dfa(y = y, num_trends = 3)
rstan::traceplot(fit3, pars = "Z")

# Illustrate how to get trends 
rotated1 = rotate_trends(fit2)
matplot(rotated1$trends)

# 
dat = read.csv("/users/eric.ward/downloads/example alaska data.csv")
dat = dat[,-1]
dat = t(dat)

fit_ak = fit_dfa(y = dat, num_trends = 3, iter=2000)
x = extract(fit_ak)$x

df = data.frame("year"=1950:2016, "mean" = apply(x[,1,], 2, mean), "low"=apply(x[,1,], 2, quantile,0.025),
  "high"=apply(x[,1,], 2, quantile,0.975), "trend"=1)
df2 = data.frame("year"=1950:2016, "mean" = apply(x[,2,], 2, mean), "low"=apply(x[,2,], 2, quantile,0.025),
  "high"=apply(x[,2,], 2, quantile,0.975), "trend"=2)
df3 = data.frame("year"=1950:2016, "mean" = apply(x[,3,], 2, mean), "low"=apply(x[,3,], 2, quantile,0.025),
  "high"=apply(x[,3,], 2, quantile,0.975), "trend"=3)
df = rbind(df, df2, df3)

ggplot(df, aes(year, mean)) + 
  geom_ribbon(aes(ymin = low, ymax=high),alpha=0.4) + 
  facet_wrap(~trend) + geom_line()

# do same model from MARSS
model.list = list(m=2, R="diagonal and unequal") 
mod3= MARSS(y, model=model.list, z.score=TRUE, form="dfa")
best.fit = mod3
H.inv = varimax(coef(best.fit, type="matrix")$Z)$rotmat

# rotate factor loadings
Z.rot = coef(best.fit, type="matrix")$Z %*% H.inv
# rotate trends
trends.rot = t(solve(H.inv) %*% best.fit$states)

marss_est = t(trends.rot)
bayes_est = t(apply(x,c(2,3),mean))


# 
ak_full = read.csv("/users/eric.ward/downloads/fate data Litzow 2-1_working.csv")
ak_full = ak_full[which(!is.na(ak_full$year)),]
meta = read.csv("/users/eric.ward/downloads/fate data Litzow 2-1_meta.csv")
meta = meta[which(meta$System != ""),]

meta = meta[which(meta$System=="GOA"),]

ak_full = ak_full[,which(names(ak_full)%in%meta$Code)]

for(i in 1:ncol(ak_full)) ak_full[,i] = as.numeric(ak_full[,i])

mod1 = fit_dfa(y = t(ak_full), num_trends = 1, iter=2000)
mod2 = fit_dfa(y = t(ak_full), num_trends = 2, iter=2000)
mod3 = fit_dfa(y = t(ak_full), num_trends = 3, iter=2000)
mod4 = fit_dfa(y = t(ak_full), num_trends = 4, iter=2000)

# calculate correlation between predicted/obs
apply(extract(mod2)$x, c(2,3), mean) %*% rotate_trends(mod2)$trends

cor(c(rotate_trends(mod2)$Z_rot %*% t(rotate_trends(mod2)$trends)), c(t(scale(ak_full))), use="pairwise.complete.obs")
cor(c(rotate_trends(mod3)$Z_rot %*% t(rotate_trends(mod3)$trends)), c(t(scale(ak_full))), use="pairwise.complete.obs")
cor(c(rotate_trends(mod4)$Z_rot %*% t(rotate_trends(mod4)$trends)), c(t(scale(ak_full))), use="pairwise.complete.obs")

trends = rotate_trends(mod2)$trends
par(mfrow = c(2,1), mai=c(0.5,0.5,0.1,0.1))
plot(1950:2016, trends[,1], type="l")
plot(1950:2016, trends[,2], type="l")

Z = rotate_trends(mod2)$Z_rot

ggplot(df, aes(year, mean)) + 
  geom_ribbon(aes(ymin = low, ymax=high),alpha=0.4) + 
  facet_wrap(~trend) + geom_line()



