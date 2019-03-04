

packageStartupMessage('Compiling model (this will take a minute...)')

dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

packageStartupMessage(paste('Writing model to:', dest))
packageStartupMessage(paste('Compiling using binary:', R.home('bin')))

model.src <- file.path(R_PACKAGE_SOURCE, 'inst', 'stan', 'corr.stan')
model.binary <- file.path(dest, 'corr_stan_model.RData')

# See: https://github.com/r-lib/pkgbuild/issues/54#issuecomment-448702834
# TODO: move stan compilation into Makevars
suppressMessages({
  model.stanc <- rstan::stanc(model.src)
  model.stanm <- rstan::stan_model(
    stanc_ret = model.stanc,
    model_name = 'corr_model'
  )
})

save('model.stanm', file = model.binary)

model.src <- file.path(R_PACKAGE_SOURCE, 'inst', 'stan', 'dfa.stan')
model.binary <- file.path(dest, 'dfa_stan_model.RData')

# See: https://github.com/r-lib/pkgbuild/issues/54#issuecomment-448702834
# TODO: move stan compilation into Makevars
suppressMessages({
  model.stanc <- rstan::stanc(model.src)
  model.stanm <- rstan::stan_model(
    stanc_ret = model.stanc,
    model_name = 'dfa_model'
  )
})

save('model.stanm', file = model.binary)

model.src <- file.path(R_PACKAGE_SOURCE, 'inst', 'stan', 'regime_1.stan')
model.binary <- file.path(dest, 'regime_1_stan_model.RData')

# See: https://github.com/r-lib/pkgbuild/issues/54#issuecomment-448702834
# TODO: move stan compilation into Makevars
suppressMessages({
  model.stanc <- rstan::stanc(model.src)
  model.stanm <- rstan::stan_model(
    stanc_ret = model.stanc,
    model_name = 'regime_1_model'
  )
})

save('model.stanm', file = model.binary)

model.src <- file.path(R_PACKAGE_SOURCE, 'inst', 'stan', 'hmm_gaussian.stan')
model.binary <- file.path(dest, 'hmm_gaussian_stan_model.RData')

# See: https://github.com/r-lib/pkgbuild/issues/54#issuecomment-448702834
# TODO: move stan compilation into Makevars
suppressMessages({
  model.stanc <- rstan::stanc(model.src)
  model.stanm <- rstan::stan_model(
    stanc_ret = model.stanc,
    model_name = 'hmm_gaussian_model'
  )
})

save('model.stanm', file = model.binary)

packageStartupMessage('------ Model successfully compiled!')
packageStartupMessage('You can ignore any compiler warnings above.')
