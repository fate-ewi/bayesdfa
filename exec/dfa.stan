data {
  int<lower=0> N; # number of data points
  int<lower=0> P; # number of time series of data
  int<lower=0> K; # number of trends
  int<lower=0> nZ; # number of unique z elements
  int<lower=0> row_indx[nZ];
  int<lower=0> col_indx[nZ];
  int<lower=0> nVariances;
  int<lower=0> varIndx[P];
  int<lower=0> nZero;
  int<lower=0> row_indx_z[nZero];
  int<lower=0> col_indx_z[nZero];
  int<lower=0> n_pos;
  real y[n_pos]; # vectorized matrix of observations
  int<lower=0> row_indx_pos[n_pos];
  int<lower=0> col_indx_pos[n_pos];
  real<lower=1> nu_fixed; // df on student-t
  int<lower=0> num_covar; # number of unique covariates
  int<lower=0> num_unique_covar; # number of covar parameters to estimate
  matrix[num_covar,N] d_covar; # inputted covariate matrix
  int covar_indexing[P,num_covar]; # index of covariates to estimate
  int estimate_nu; # Estimate degrees of freedom?
}
parameters {
  matrix[K,N] x; #vector[N] x[P]; # random walk-trends
  vector[nZ] z; # estimated loadings in vec form
  real<lower=0> sigma[nVariances];
  // real<lower=2> nu;
  real<lower=2> nu[estimate_nu]; // df on student-t
}
transformed parameters {
  matrix[P,N] pred; #vector[P] pred[N];
  matrix[P,K] Z;
  for(i in 1:nZ) {
    Z[row_indx[i],col_indx[i]] = z[i];
  }
  # fill in zero elements
  if(nZero > 2) {
    for(i in 1:(nZero-2)) {
      Z[row_indx_z[i],col_indx_z[i]] = 0;
    }
  }

  for(k in 1:K) {
    Z[k,k] = 1;// add constraint for Z diagonal
  }
  # N is sample size, P = time series, K = number trends
  # [PxN] = [PxK] * [KxN]
  pred = Z * x;
}
model {
  # initial state for each trend
  for(k in 1:K) {
    x[k,1] ~ normal(0, 1);
    for(t in 2:N) {
      if (estimate_nu == 1) {
        x[k,t] ~ student_t(nu[1], x[k,t-1], 1); # random walk
      } else {
        x[k,t] ~ student_t(nu_fixed, x[k,t-1], 1); # random walk
      }
    }
  }
  if (estimate_nu == 1) {
    nu[1] ~ gamma(2, 0.1);
    // nu[1] ~ exponential(0.01);
  }

  # prior on loadings
  z ~ normal(0, 1);

  # observation variance
  sigma ~ student_t(3, 0, 2);

  # likelihood
  for(i in 1:n_pos) {
    y[i] ~ normal(pred[row_indx_pos[i], col_indx_pos[i]], sigma[varIndx[row_indx_pos[i]]]);
  }

}
generated quantities {
  vector[n_pos] log_lik;
   #regresssion example in loo() package
  for (n in 1:n_pos) {
    log_lik[n] = normal_lpdf(y[n] | pred[row_indx_pos[n], col_indx_pos[n]], sigma[varIndx[row_indx_pos[n]]]);
  }
}
