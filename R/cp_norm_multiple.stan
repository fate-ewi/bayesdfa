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
  real m;
  real<lower=0> sigma;
}
transformed parameters {
  matrix[T, T] lp;
  lp = rep_matrix(log_unif, T, T);
  for (s1 in 1:T)
    for (s2 in 1:T)
      for (t in 1:T)
        lp[s1,s2] = lp[s1,s2] +
          + normal_lpdf(D[t] | t < s1 ? e : (t < s2 ? m: l), sigma);
}
model {
  e ~ normal(0, 2);
  m ~ normal(0, 2);
  l ~ normal(0, 2);
  sigma ~ student_t(3, 0, 1);
  target += log_sum_exp(to_vector(lp));
}
generated quantities {
  int<lower=1,upper=T> s;
  vector[T] expectations;
  s = categorical_logit_rng(to_vector(lp));
  expectations = softmax(to_vector(lp));
}
