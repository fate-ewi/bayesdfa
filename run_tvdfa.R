library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")
data(harborSealWA)
y <- t(harborSealWA[, c("SJF", "SJI", "EBays", "PSnd", "HC")])

m <- fit_dfa(y = y, num_trends = 1, iter = 1000, model = "tvdfa.stan")
m_ntv <- fit_dfa(y = y, num_trends = 1, iter = 2000, model = "dfa.stan")

# fit1
b <- broom::tidyMCMC(m, rhat = TRUE, ess = TRUE, conf.int = TRUE, conf.method = "quantile")
b <- filter(b, !is.na(rhat)) # fixed at 0
max(b$rhat)
min(b$ess)

library(tidyverse)
z <- filter(b, grepl("Z\\[", term))
z <- mutate(z, 
  time = as.numeric(gsub("Z\\[([0-9]+),[0-9]+,[0-9]+]", "\\1", z$term)),
  ts = gsub("Z\\[([0-9])+,([0-9]+),[0-9]+]", "\\2", z$term),
  trend = gsub("Z\\[([0-9])+,([0-9])+,([0-9]+)]", "\\3", z$term)
)
ggplot(z, aes(time, estimate, colour = as.factor(ts))) + geom_line() +
  theme_light() +
  facet_wrap(~trend) #+
  # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(ts)), alpha = 0.1, lwd = 0)

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
b2 <- filter(b2, !is.na(rhat)) # fixed at 0
max(b2$rhat)
min(b2$ess)

z <- filter(b2, grepl("Z\\[", term))
z <- mutate(z, 
  ts = as.numeric(gsub("Z\\[([0-9])+,([0-9]+)]", "\\1", z$term)),
  trend = as.numeric(gsub("Z\\[([0-9])+,([0-9]+)]", "\\2", z$term))
)
# ggplot(z, aes(time, estimate, colour = as.factor(ts))) + geom_bar() +
  # theme_light() +
  # facet_wrap(~trend)


x <- filter(b2, grepl("x\\[", term))
x <- mutate(x, 
  time = as.numeric(gsub("x\\[([0-9]+),([0-9]+)]", "\\2", x$term)),
  trend = as.numeric(gsub("x\\[([0-9]+),([0-9])+]", "\\1", x$term))
)
ggplot(x, aes(time, estimate, colour = as.factor(trend))) + geom_line() +
  theme_light() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = as.factor(trend)), alpha = 0.1, lwd = 0)
