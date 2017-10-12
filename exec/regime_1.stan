data {
  int<lower=0> N; // number of data points
  vector[N] y; // data, no missing values
}
parameters {
  real u;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] pred;
  for(i in 1:N) {
    pred[i] = u; // default
  }
}
model {
  // priors
  sigma ~ student_t(3, 0, 2);
  u ~ normal(0, 1);
  y ~ normal(pred, sigma);
}
generated quantities {
  vector[N] log_lik;
   //regresssion example in loo() package
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | pred[n], sigma);
  }
}
