data {
  int<lower=0> N; # number of data points
  int<lower=0> P; # number of time series of data
  int<lower=0> K; # number of trends
  int<lower=0> nZ; # number of unique z elements
  vector[N] y[P]; # matrix of observations
  int<lower=0> row_indx[nZ];
  int<lower=0> col_indx[nZ];
  int<lower=0> nVariances;
  int<lower=0> varIndx[P];
  int<lower=0> nZero;  
  int<lower=0> row_indx_z[nZero];
  int<lower=0> col_indx_z[nZero];  
}
parameters {
  matrix[K,N] x; #vector[N] x[P]; # random walk-trends
  vector[nZ] z; # estimated loadings in vec form
  real<lower=0> sigma[nVariances];
}
transformed parameters {
  matrix[P,N] pred; #vector[P] pred[N];
  matrix[P,K] Z;
  for(i in 1:nZ) {
    Z[row_indx[i],col_indx[i]] = z[i];
  }
  # fill in zero elements
  #for(i in 1:nZero) {
  #  Z[row_indx_z[i],col_indx_z[i]] = 0;
  #}
  Z[1,2] = 0;
  
  # N is sample size, P = time series, K = number trends
  # [PxN] = [PxK] * [KxN]
  pred = Z * x;
}
model {
  # initial state for each trend
  for(k in 1:K) {
    x[k,1] ~ normal(0,3);
    for(t in 2:N) {
      x[k,t] ~ normal(x[k,t-1],1); # random walk
    }
  }
  # prior on loadings
  for(i in 1:nZ) {z[i] ~ normal(0,1);}
  
  # observation variance
  for(i in 1:nVariances) {
  sigma[i] ~ cauchy(0, 5);
  }
  
  # likelihood
  for(p in 1:P) {
  y[p,] ~ normal(pred[p,], sigma[varIndx[p]]);
  }


}

