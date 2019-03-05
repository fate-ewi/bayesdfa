data {
  int<lower=1> T;                   // number of observations (length)
  int<lower=1> K;                   // number of hidden states
  real x_t[T];                      // observations
  int<lower=0> est_sigma;           // flag, whether to estimate sigma (1) or use values passed in (0)
  real sigma_t[T];               // estimated sigma for each observation
}
parameters {
  real mu_k;                  // observation means
  real<lower=0> sigma_k;         // observation standard deviations, optionally estimated if est_sigma == 1. Can the quantity K * est_sigma be used to dimension sigma_k?
}
transformed parameters {
  real sigmas[T];
  if(est_sigma == 1) {
    for(i in 1:T) sigmas[i] = sigma_k;
  } else {
    for(i in 1:T) sigmas[i] = sigma_t[i];
  }
}
model {
  mu_k ~ student_t(3, 0, 3);
  sigma_k ~ student_t(3, 0, 1);

  x_t ~ normal(mu_k, sigmas);
}
generated quantities {
  vector[T] log_lik;
   //regresssion example in loo() package
  for (n in 1:T) {
    log_lik[n] = normal_lpdf(x_t[n] | mu_k, sigmas[n]);
  }
}
