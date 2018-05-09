library(bayesdfa)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# These are all the simulation parameters
mcmc_chains = 3
mcmc_iter = 4000
mcmc_warm = 3000
sim_years = 20
sim_ts = 6

set.seed(123)

iter = 100
df = data.frame(seeds = sample.int(.Machine$integer.max, iter),
  conv = NA,
  trnd_sim = sample(1:4, size=iter, replace=T),
  num_trends = sample(1:4, size=iter,replace=T),
  sigma = sample(c(0.1,1),size=iter, replace=T),
  pars_not_conv = NA,
  max_rhat = NA
  )

# Run models on simulated data 
for(ii in 1:iter) {

set.seed(df$seeds[ii])
# Simulate DFA data, varying trends, and obs error
y = bayesdfa::sim_dfa(num_trends = df$trnd_sim[ii], 
  num_years = sim_years, num_ts=sim_ts, sigma=df$sigma[ii])$y_sim

# Fit a different model to the data. Data generating
# models almost always work, but we're more interested
# in cases of poor convergence
fit_t = fit_dfa(y = y, 
  num_trends = df$num_trends[ii], 
  seed = df$seeds[ii],
  iter=mcmc_iter, 
  warmup = mcmc_warm, 
  chains=mcmc_chains,
  varIndx = rep(1, nrow(y)))
# Invert the chains
invert = invert_chains(fit_t$model, trends=df$num_trends[ii])
# Evaluate convergence
df$conv[ii] = is_converged(invert)
# Calculate maximum Rhat and number of parameters not converged
rhats = as.data.frame(invert$monitor)[,"Rhat"]
df$max_rhat[ii] = max(rhats,na.rm=T)
df$pars_not_conv[ii] = length(which(rhats > 1.05))
}

saveRDS(df, "test_dfa.rds")
