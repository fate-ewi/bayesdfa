data {
  int<lower=1> T;
  real D[T];
}
transformed data {
  real log_unif;
  log_unif = -log(T);
}
parameters {
  real e;
  real l;
  real<lower=0> sigma;
}
transformed parameters {
  vector[T] lp;
  lp = rep_vector(log_unif, T);
  for (s in 1:T)
    for (t in 1:T)
      lp[s] = lp[s] + normal_lpdf(D[t] | t < s ? e : l, sigma);
}
model {
  e ~ normal(0, 2);
  l ~ normal(0, 2);
  sigma ~ student_t(3, 0, 1);
  target += log_sum_exp(lp);
}
generated quantities {
  int<lower=1,upper=T> s;
  vector[T] expectations;
  s = categorical_logit_rng(lp);
  expectations = softmax(lp);
}
