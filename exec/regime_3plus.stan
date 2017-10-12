data {
  int<lower=0> N; // number of data points
  int<lower=0> n_regime; // dimension of change points (this is number of regimes)
  vector[N] y; // data, no missing values
  vector[(n_regime-1)] ones;
}
parameters {
  vector[n_regime] u;
  simplex[(n_regime-1)] theta;
  real<lower=0> sigma;
}
transformed parameters {
  vector[N] pred;
  vector[(n_regime-1)] x;// estimated locations of change points
  vector[(n_regime-1)] sumtheta;

  sumtheta[1] = theta[1];
  for(i in 2:(n_regime-1)) {
    sumtheta[i] = sumtheta[i-1] + theta[i]; // cumulative, sum to 1
  }
  for(i in 1:n_regime-1) {
    x[i] = 1 + (N-1) * sumtheta[i];
  }

  // assign each observation to regime/
  for(i in 1:N) {
    pred[i] = u[n_regime]; // default
    if(i < x[1]) pred[i] = u[1];
    for(j in 2:(n_regime-1)) {
      if(i > x[j-1] && i < x[j]) pred[i] = u[j];
    }
  }
}
model {
  // priors
  theta ~ dirichlet(ones);
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
