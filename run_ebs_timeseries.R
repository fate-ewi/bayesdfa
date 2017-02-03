library(readr)
meta <- read_csv("C:/Users/great/Desktop/EWI_Nonlinearity/Alaska/fate data Litzow 2-1_meta.csv")
working <- read_csv("C:/Users/great/Desktop/EWI_Nonlinearity/Alaska/fate data Litzow 2-1_working.csv")

Data_ori<-working[1:66,]
Meta<-meta[1:79,]


j=1
stocks=rep(NA,79)
for (i in 1:79) {
  if ((grepl("CI",Meta$Name[i])==FALSE)&&(grepl("sample size",Meta$Name[i])==FALSE)) {
    if (Meta$System[i]=="EBS"){
      stocks[j]=Meta$Code[i]
      j=j+1
    }
  }
}

EBS_ind=match(stocks, colnames(Data_ori))
ebs_full<-Data_ori[,EBS_ind[1:length(which(!is.na(stocks)))]]
rm(list=setdiff(ls(), "ebs_full"))
ebs<-as.matrix(ebs_full)
ebs_full<-matrix(NA,dim(ebs)[1],dim(ebs)[2])
for(i in 1:ncol(ebs_full)) ebs_full[,i] = as.numeric(ebs[,i])

##########run stan
library(rstan)
library(MARSS)
library(loo)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
source("find_dfa_trends.r")

df = find_dfa_trends(y = t(ebs_full), iter=200)


mod1 = fit_dfa(y = t(ebs_full), num_trends = 1, iter=2000, varIndx = seq(1,ncol(ebs_full)))
rt = rotate_trends(mod1)

par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(1950:2015, apply(rt$trends[,1,],2,mean), type="l", ylim=c(-4,4), lwd=3)
lines(1950:2015, apply(rt$trends[,1,],2,quantile,0.025,col="grey"))
lines(1950:2015, apply(rt$trends[,1,],2,quantile,0.975,col="grey"))


par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(1950:2014, diff(apply(rt$trends[,1,],2,mean)), type="l", lwd=3, ylab="anomaly")
#plot(1950:2015, diff(apply(rt$trends[,2,],2,mean)), type="l", lwd=3, ylab="anomaly")
#plot(1950:2015, diff(apply(rt$trends[,3,],2,mean)), type="l", lwd=3, ylab="anomaly")
#xx = cbind(diff(apply(rt$trends[,1,],2,mean)), diff(apply(rt$trends[,2,],2,mean)), diff(apply(rt$trends[,3,],2,mean)))
#plot(1950:2015, log(mvtnorm::dmvnorm(xx, mean = rep(0,3))), type="l", ylab="log(dmvnorm)")

par(mfrow = c(2,2), mai=c(0.5,0.7,0.1,0.1))
plot(rollapply(zoo(apply(rt$trends[,1,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")
#plot(rollapply(zoo(apply(rt$trends[,2,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")
#plot(rollapply(zoo(apply(rt$trends[,3,],2,sd)), 3, sd), type="l", lwd=3, ylab="sd")

#Z = rotate_trends(mod1)$Z_rot

#Zmean = cbind(colnames(ebs_full), round(apply(Z, c(2,3), mean), 2))

#ggplot(df, aes(year, mean)) + 
#geom_ribbon(aes(ymin = low, ymax=high),alpha=0.4) + 
#facet_wrap(~trend) + geom_line()


