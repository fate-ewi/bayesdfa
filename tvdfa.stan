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
}
parameters {
  matrix[K,N] x; #vector[N] x[P]; # random walk-trends
  vector[nZ] z[N]; # estimated loadings in vec/matrix form
  real<lower=0> sigma[nVariances];
  real<lower=0> tau[P,K];
}
transformed parameters {
  matrix[P,N] pred;
  matrix[P,K] Z[N];
  
  // print("Z=", Z);
  for(n in 1:N) {
    for(i in 1:nZ) {
      Z[n, row_indx[i], col_indx[i]] = z[n][i];
    }
  }

  
  # fill in zero elements
  if(nZero > 2) {
    for(n in 1:N) {
      for(i in 1:(nZero-2)) {
        Z[n, row_indx_z[i],col_indx_z[i]] = 0;
      }
    }
  }
  
  # N is sample size, P = time series, K = number trends
  # [PxN] = [PxK] * [KxN]
  for(n in 1:N) {
    for(p in 1:P) {
        pred[p,n] = Z[n,p,1] * x[1,n];
    }
  }
  if (K > 1) {
    for(n in 1:N) {
      for(p in 1:P) {
        for(k in 2:K) {
          pred[p,n] = pred[p,n] + Z[n,p,k] * x[k,n];
        }
      }
    }
  }
}
model {
  # initial state for each factor trend
  for(k in 1:K) {
    x[k,1] ~ normal(0, 1);
    for(t in 2:N) {
      x[k,t] ~ normal(x[k,t-1], 1); # random walk
    }
  }

  # initial state for each load trend
  for(k in 1:K) {
    for(p in 1:P) {
      Z[1,p,k] ~ normal(0, 1);
    }
  }
  for(n in 2:N) {
    for(k in 1:K) {
      for(p in 1:P) {
        Z[n,p,k] ~ normal(Z[n-1,p,k], tau[p,k]); # random walk
      }
    }
  }

  for(k in 1:K) {
    tau[k] ~ student_t(3, 0, 2);
  }
  
  # prior on loadings
  for(n in 1:N) {
    z[n] ~ normal(0, 1);
  }
  
  # observation variance
  for(i in 1:nVariances) {
    sigma[i] ~ student_t(3, 0, 2);
  }
  
  # likelihood
  for(i in 1:n_pos) {
    y[i] ~ normal(pred[row_indx_pos[i], col_indx_pos[i]], sigma[varIndx[row_indx_pos[i]]]);
  }

}
generated quantities {
  vector[n_pos] log_lik;
  # regresssion example in loo() package 
  for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[row_indx_pos[n], col_indx_pos[n]], sigma[varIndx[row_indx_pos[n]]]);
}
