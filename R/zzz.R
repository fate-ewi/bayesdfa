## Copyright (c) 2017-present, Facebook, Inc.
## All rights reserved.

## This source code is licensed under the BSD-style license found in the
## LICENSE file in the root directory of this source tree. An additional grant
## of patent rights can be found in the PATENTS file in the same directory.

.onLoad <- function(libname, pkgname) {
  .dfa.stan.model <- get_stan_model(model_name="dfa")
  assign(
    ".dfa.stan.model",
    .dfa.stan.model,
    envir=parent.env(environment())
  )

  .corr.stan.model <- get_stan_model(model_name="corr")
  assign(
    ".corr.stan.model",
    .corr.stan.model,
    envir=parent.env(environment())
  )

  .hmm.stan.model <- get_stan_model(model_name="hmm_gaussian")
  assign(
    ".hmm.stan.model",
    .hmm.stan.model,
    envir=parent.env(environment())
  )

  .regime.stan.model <- get_stan_model(model_name="regime_1")
  assign(
    ".regime.stan.model",
    .regime.stan.model,
    envir=parent.env(environment())
  )
}
