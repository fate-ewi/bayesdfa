data {
  int<lower=0> N; // number of data points
  int<lower=0> P; // number of time series of data
  int<lower=0> K; // number of trends
  int<lower=0> nZ; // number of unique z elements
  int<lower=0> row_indx[nZ];
  int<lower=0> col_indx[nZ];
  int<lower=0> nVariances;
  int<lower=0> varIndx[P];
  int<lower=0> nZero;
  int<lower=0> row_indx_z[nZero];
  int<lower=0> col_indx_z[nZero];
  int<lower=0> n_pos; // number of non-missing observations
  int<lower=0> row_indx_pos[n_pos]; // row indices of non-missing obs
  int<lower=0> col_indx_pos[n_pos]; // col indices of non-missing obs
  real y[n_pos]; // vectorized matrix of observations
  int<lower=0> n_na; // number of missing observations
  int<lower=0> row_indx_na[n_na]; // row indices of missing obs
  int<lower=0> col_indx_na[n_na]; // col indices of missing obs
  real<lower=1> nu_fixed; // df on student-t
  int<lower=0> num_covar; // number of unique covariates
  int<lower=0> num_unique_covar; // number of covar parameters to estimate
  matrix[num_covar,N] d_covar; // inputted covariate matrix
  int covar_indexing[P,num_covar]; // index of covariates to estimate
  int estimate_nu; // Estimate degrees of freedom?
  int use_normal; // flag, for large values of nu > 100, use normal instead
  int est_cor; // whether to estimate correlation in obs error (=1) or not (=0)
  int est_phi; // whether to estimate autocorrelation in trends (=1) or not (= 0)
}
transformed data {
  int n_pcor; // dimension for cov matrix
  int n_loglik; // dimension for loglik calculation

  if(est_cor == 0) {
    n_loglik = P;
  } else {
    n_loglik = N;
  }

  if(est_cor == 0) {
    n_pcor = P;
    if(nVariances < 2) {
      n_pcor = 2;
    }
  } else {
    n_pcor = P;
  }
}
parameters {
  matrix[K,N] x; //vector[N] x[P]; // random walk-trends
  vector[nZ] z; // estimated loadings in vec form
  vector<lower=0>[K] zpos; // constrained positive values
  real<lower=0> sigma[nVariances];
  real<lower=2> nu[estimate_nu]; // df on student-t
  real ymiss[n_na];
  real<lower=-1,upper=1> phi[est_phi*K];
  cholesky_factor_corr[n_pcor] Lcorr;
}
transformed parameters {
  matrix[P,N] pred; //vector[P] pred[N];
  matrix[P,K] Z;
  //vector[N] yall[P]; // combined vectors of missing and non-missing values
  matrix[P,N] yall;
  vector[P] sigma_vec;
  vector[P] phi_vec;

  if(est_phi == 1) {
    for(p in 1:P) {phi_vec[p] = phi[p];}
  } else {
    for(p in 1:P) {phi_vec[p] = 1;}
  }

  for(p in 1:P) {
    sigma_vec[p] = sigma[varIndx[p]]; // convert estimated sigmas to vec form
  }

  // Fill yall with non-missing values
  for(i in 1:n_pos) {
    yall[row_indx_pos[i], col_indx_pos[i]] = y[i];
  }
  // Include missing observations
  if(n_na > 0) {
    for(i in 1:n_na) {
      yall[row_indx_na[i], col_indx_na[i]] = ymiss[i];
    }
  }

  for(i in 1:nZ) {
    Z[row_indx[i],col_indx[i]] = z[i]; // convert z to from vec to matrix
  }
  // fill in zero elements
  if(nZero > 2) {
    for(i in 1:(nZero-2)) {
      Z[row_indx_z[i],col_indx_z[i]] = 0;
    }
  }

  for(k in 1:K) {
    Z[k,k] = zpos[k];// add constraint for Z diagonal
  }
  // N is sample size, P = time series, K = number trends
  // [PxN] = [PxK] * [KxN]
  pred = Z * x;
}
model {
  // initial state for each trend
  for(k in 1:K) {
    x[k,1] ~ cauchy(0,3);//normal(0, 1);
    if(use_normal == 0) {
      for(t in 2:N) {
        if (estimate_nu == 1) {
          x[k,t] ~ student_t(nu[1], phi_vec[k]*x[k,t-1], 1); // random walk
        } else {
          x[k,t] ~ student_t(nu_fixed, phi_vec[k]*x[k,t-1], 1); // random walk
        }
      }
    } else {
      for(t in 2:N) {
        x[k,t] ~ normal(phi_vec[k]*x[k,t-1], 1);
      }
    }

  }
  if (estimate_nu == 1) {
    nu[1] ~ gamma(2, 0.1); // df parameter for t-distribution
  }

  if(est_phi == 1) {
    for(p in 1:P) {
      // uniform prior on AR prior if included
      phi[p] ~ uniform(-1,1);
    }
  }
  // prior on loadings
  z ~ normal(0, 1);
  zpos ~ normal(0, 1);

  // observation variance
  sigma ~ student_t(3, 0, 2);
  if(est_cor == 1) {
    Lcorr ~ lkj_corr_cholesky(1);
  }
  // likelihood for independent
  if(est_cor == 0) {
    for(i in 1:P){
      target += normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
    }
  } else {
    // need to loop over time slices / columns - each ~ MVN
    for(i in 1:N) {
      target += multi_normal_cholesky_lpdf(col(yall,i) | col(pred,i), diag_pre_multiply(sigma_vec, Lcorr));
    }
  }

}
generated quantities {
  vector[n_loglik] log_lik;
  matrix[n_pcor, n_pcor] Omega;
  matrix[n_pcor, n_pcor] Sigma;
  if(est_cor == 1) {
  Omega = multiply_lower_tri_self_transpose(Lcorr);
  Sigma = quad_form_diag(Omega, sigma_vec);
  }

  //calculate looic based on regresssion example in loo() package
  if(est_cor == 0) {
    //for (n in 1:n_pos) {
    //  log_lik[n] = normal_lpdf(y[n] | pred[row_indx_pos[n], col_indx_pos[n]], sigma[varIndx[row_indx_pos[n]]]);
    //}
    for(i in 1:P) {
      log_lik[i] = normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
    }
  } else {
    for(i in 1:N) {
      log_lik[i] = multi_normal_cholesky_lpdf(col(yall,i) | col(pred,i), diag_pre_multiply(sigma_vec, Lcorr));
    }
  }
}
