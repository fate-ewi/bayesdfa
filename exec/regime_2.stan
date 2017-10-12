data {
  int<lower=0> N; // number of data points
  vector[N] y; // data, no missing values
}
parameters {
  vector[2] u;
  real<lower=0,upper=1> theta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] pred;
  real x;// estimated locations of change points

  x = 1 + (N-1) * theta;

  // assign each observation to regime/
  for(i in 1:N) {
    pred[i] = u[2]; // default
    if(i < x) pred[i] = u[1];
  }
}
model {
  // priors
  theta ~ beta(1,1);
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
