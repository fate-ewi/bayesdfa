data {
  int<lower=0> N;
  vector[N] y;
  vector[N] x;
}
parameters {
  real<lower=-1,upper=1> beta;
  real<lower=0> sigma;
}
model {
  // priors
  sigma ~ student_t(3, 0, 2);
  y ~ normal(beta * x, sigma);
}
