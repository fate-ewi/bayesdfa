library(rstan)
library(MARSS)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("fit_dfa.r")
source("rotate_trends.r")


# load the data (there are 3 datasets contained here)
 data(lakeWAplankton)
 # we want lakeWAplanktonTrans, which has been transformed
 # so the 0s are replaced with NAs and the data z-scored
 dat = lakeWAplanktonTrans
 # use only the 10 years from 1980-1989
 plankdat = dat[dat[,"Year"]>=1980 & dat[,"Year"]<1990,]
 # create vector of phytoplankton group names
 phytoplankton = c("Cryptomonas", "Diatoms", "Greens",
                   "Unicells", "Other.algae")
 # get only the phytoplankton
 dat.spp.1980 = plankdat[,phytoplankton]

# fit the three trend model with 3 covariates
y = t(dat.spp.1980)
covar = matrix(runif(3*ncol(y)), nrow=3)
fit_dfa(y = y, covar=covar, num_trends=3)

# fit the three trend model with 37 covariates
y = t(dat.spp.1980)
covar = matrix(runif(37*ncol(y)), nrow=37)
covar_index = matrix(1, nrow(y), nrow(covar))
# let's make species 1, 3 have the same effects, 2/4 have the same effects and 5 be different
covar_index[1,] = 1:37
covar_index[3,] = 1:37
covar_index[2,] = 38:74
covar_index[4,] = 38:74
covar_index[5,] = 75:111
fit_dfa(y = y, covar=covar, covar_index = covar_index, num_trends=3)


