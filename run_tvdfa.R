library(rstan)
library(MARSS)
library(tidyverse)
rstan_options(auto_write = TRUE)

source("fit_dfa.r")
source("rotate_trends.r")
data(harborSealWA)
d <- t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd")])

d <- as.data.frame(ivesDataLP[,c("LargePhyto", "SmallPhyto")])
d$time <- 1:nrow(d)
gather(d, size, value, - time) %>%
  ggplot(aes(time, value, colour = size)) + geom_line()
d <- t(ivesDataLP[,c("LargePhyto", "SmallPhyto")])

m <- fit_dfa(y = log(d), num_trends = 1, iter = 3000, model = "tvdfa.stan", varIndx = 1:nrow(d))
m_ntv <- fit_dfa(y = log(d), num_trends = 1, iter = 4000, model = "dfa.stan", varIndx = 1:nrow(d),
  nu = 1e6)
m_ntv_nu1 <- fit_dfa(y = log(d), num_trends = 1, iter = 4000, model = "dfa.stan", varIndx = 1:nrow(d),
  nu = 3)

library(loo)
loo(extract_log_lik(m))$looic
loo(extract_log_lik(m_ntv))$looic

# fit1
b <- broom::tidyMCMC(m, rhat = TRUE, ess = TRUE, conf.int = TRUE, conf.method = "quantile")
b <- filter(b, !is.na(rhat)) # fixed at 0
max(b$rhat)
min(b$ess)
filter(b, grepl("tau", term))
filter(b, grepl("sigma", term))

z <- filter(b, grepl("Z\\[", term))
z <- mutate(z, 
  time = as.numeric(gsub("Z\\[([0-9]+),[0-9]+,[0-9]+]", "\\1", z$term)),
  ts = gsub("Z\\[([0-9])+,([0-9]+),[0-9]+]", "\\2", z$term),
  trend = gsub("Z\\[([0-9])+,([0-9])+,([0-9]+)]", "\\3", z$term)
)
ggplot(z, aes(time, estimate, colour = as.factor(ts))) + geom_line() +
  theme_light() +
  facet_wrap(~trend) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(ts)), alpha = 0.1, lwd = 0)

x <- filter(b, grepl("x\\[", term))
x <- mutate(x, 
  time = as.numeric(gsub("x\\[([0-9]+),([0-9]+)]", "\\2", x$term)),
  trend = as.numeric(gsub("x\\[([0-9]+),([0-9])+]", "\\1", x$term))
)
ggplot(x, aes(time, estimate, colour = as.factor(trend))) + geom_line() +
  theme_light() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(trend)), alpha = 0.1, lwd = 0)

#####
b2 <- broom::tidyMCMC(m_ntv, rhat = TRUE, ess = TRUE, conf.int = TRUE, conf.method = "quantile")
b2_t <- broom::tidyMCMC(m_ntv_nu1, rhat = TRUE, ess = TRUE, conf.int = TRUE, conf.method = "quantile")
b2 <- filter(b2, !is.na(rhat)) # fixed at 0
max(b2$rhat)
min(b2$ess)

z2 <- filter(b2, grepl("Z\\[", term))
z2 <- mutate(z2, 
  ts = gsub("Z\\[([0-9])+,([0-9]+)]", "\\1", z2$term),
  trend = gsub("Z\\[([0-9])+,([0-9]+)]", "\\2", z2$term)
)
ggplot(z2, aes(ts, estimate)) + geom_bar(stat = "identity") +
  theme_light() +
  facet_wrap(~trend)

x2 <- filter(b2, grepl("x\\[", term))
x2 <- mutate(x2, 
  time = as.numeric(gsub("x\\[([0-9]+),([0-9]+)]", "\\2", x2$term)),
  trend = as.numeric(gsub("x\\[([0-9]+),([0-9])+]", "\\1", x2$term))
)

x2_t <- filter(b2_t, grepl("x\\[", term))
x2_t <- mutate(x2_t,
  time = as.numeric(gsub("x\\[([0-9]+),([0-9]+)]", "\\2", x2_t$term)),
  trend = as.numeric(gsub("x\\[([0-9]+),([0-9])+]", "\\1", x2_t$term))
)

if (sign(x2_t$estimate[1]) != sign(x2$estimate[1])) x2_t$estimate <- x2_t$estimate * -1

ggplot(x2, aes(time, estimate, colour = as.factor(trend))) + geom_line() +
  theme_light() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(trend)), alpha = 0.1, lwd = 0) +
  geom_line(data = x2_t, aes(time, estimate), lty = 2, col = "black", inherit.aes = FALSE)
