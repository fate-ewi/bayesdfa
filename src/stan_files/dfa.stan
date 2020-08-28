functions {

  // this function subsets a matrix by dropping the row/column labeled 'drop'. P represents dimensions
  matrix subset(matrix x, int drop, int P) {
    // count number of rows in result

    // assign rows in result
    {
      matrix[P-1,P-1] result;

      int rowindx;
      int colindx;
      rowindx = 0;
      for (i in 1:P) {
        if (i != drop) {
          rowindx = rowindx + 1;
          colindx = 0;
          for (j in 1:P) {
            if (j != drop) {
              colindx = colindx + 1;
              result[rowindx, colindx] = x[i, j];
            }
          } // end j loop
        } // end i!= drop
      } // end i loop

      return result;
    }
  }

  matrix subsetvec(matrix x, int drop, int P) {
    // assign rows in result
    {
      matrix[P-1,1] result;

      int rowindx;
      rowindx = 0;
      for (i in 1:P) {
        if (i != drop) {
          rowindx = rowindx + 1;
          result[rowindx,1] = x[i, drop];
        } // end i!= drop
      } // end i loop

      return result;
    }
  }

  matrix subsetvec2(vector x, int drop, int P) {
    // assign rows in result
    {
      matrix[P-1,1] result;

      int rowindx;
      rowindx = 0;
      for (i in 1:P) {
        if (i != drop) {
          rowindx = rowindx + 1;
          result[rowindx,1] = x[i];
        } // end i!= drop
      } // end i loop

      return result;
    }
  }
}
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
  int estimate_nu; // Estimate degrees of freedom?
  int use_normal; // flag, for large values of nu > 100, use normal instead
  int est_cor; // whether to estimate correlation in obs error (=1) or not (=0)
  int est_phi; // whether to estimate autocorrelation in trends (=1) or not (= 0)
  int est_theta; // whether to estimate moving-average in trends (=1) or not (= 0
  int<lower=0> num_obs_covar; // number of unique observation covariates, dimension of matrix
  int<lower=0> n_obs_covar; // number of unique covariates included
  int obs_covar_index[num_obs_covar,3];// indexed by time, trend, covariate #, covariate value. +1 because of indexing issues
  real obs_covar_value[num_obs_covar];
  int<lower=0> num_pro_covar; // number of unique process covariates, dimension of matrix
  int<lower=0> n_pro_covar; // number of unique process covariates included
  int pro_covar_index[num_pro_covar,3] ;// indexed by time, trend, covariate #, covariate value. +1 because of indexing issues
  real pro_covar_value[num_pro_covar];
  real z_bound[2];
  int<lower=0> long_format; // data shape, 0 == wide (default), 1 = long with potential for multiple observations
  int<lower=0> proportional_model;
  int<lower=0> est_sigma_process; // optional, 0 == not estimate sigma_pro (default), 1 == estimate
  int<lower=0> n_sigma_process; // single value, or equal number of tre
}
transformed data {
  int n_pcor; // dimension for cov matrix
  int n_loglik; // dimension for loglik calculation
  vector[K] zeros;
  vector[K*proportional_model] alpha_vec;

  for(k in 1:K) {
    zeros[k] = 0; // used in MVT / MVN below
  }

  // this is number of points of log-likelihood, depends if model is MVN or not and data is in wide/long format
  n_loglik = n_pos;
  if(long_format==0) {
    if(est_cor == 0) {
       n_loglik = P * N;
    } else {
      n_loglik = N;
    }
  }

  if(est_cor == 0) {
    // MVN correlation matrix not estimated
    n_pcor = P;
    if(nVariances < 2) {
      // minimum bound, just for Stan types
      n_pcor = 2;
    }
  } else {
    // MVN correlation matrix is estimated
    n_pcor = P;
  }

  // for compositional model
  if(proportional_model==1) {
    for(k in 1:K) alpha_vec[k] = 1;
  }
}
parameters {
  matrix[K,N-1] devs; // random deviations of trends
  vector[K] x0; // initial state
  vector<lower=0>[K*(1-proportional_model)] psi; // expansion parameters
  vector<lower=z_bound[1],upper=z_bound[2]>[nZ*(1-proportional_model)] z; // estimated loadings in vec form
  vector[K*(1-proportional_model)] zpos; // constrained positive values
  simplex[K] p_z[P*proportional_model]; // alternative for proportional Z
  matrix[n_obs_covar, P] b_obs; // coefficients on observation model
  matrix[n_pro_covar, K] b_pro; // coefficients on process model
  real<lower=0> sigma[nVariances];
  real<lower=2> nu[estimate_nu]; // df on student-t
  real ymiss[n_na];
  real<lower=-1,upper=1> phi[est_phi*K];
  real<lower=-1,upper=1> theta[est_theta*K];
  cholesky_factor_corr[n_pcor] Lcorr;
  real<lower=0> sigma_process[est_sigma_process * n_sigma_process];
}
transformed parameters {
  matrix[P,N] pred; //vector[P] pred[N];
  matrix[P,K] Z;
  //vector[N] yall[P]; // combined vectors of missing and non-missing values
  matrix[P,N] yall;
  vector[P] sigma_vec;
  vector[K] phi_vec; // for AR(1) part
  vector[K] theta_vec; // for MA(1) part
  matrix[K,N] x; //vector[N] x[P]; // random walk-trends
  vector[K] indicator; // indicates whether diagonal is neg or pos
  vector[K] psi_root; // derived sqrt(expansion parameter psi)
  matrix[n_pcor*long_format*est_cor, n_pcor*long_format*est_cor] Sigma_derived;// temporary for calculations for MVN model
  matrix[(n_pcor-1)*long_format*est_cor, (n_pcor-1)*long_format*est_cor] Sigma_temp;// temporary for calculations for MVN model
  matrix[n_pcor-1,1] sigma12_vec;// temporary for calculations for MVN model
  matrix[P*long_format*est_cor, N*long_format*est_cor] temp_sums;// temporary for calculations for MVN model
  matrix[P*long_format*est_cor, N*long_format*est_cor] temp_counts;// temporary for calculations for MVN model
  vector[P*long_format*est_cor] cond_sigma_vec;// temporary for calculations for MVN model
  vector[P*long_format*est_cor] cond_mean_vec;// temporary for calculations for MVN model
  real sigma11;// temporary for calculations for MVN model

  // phi is the ar(1) parameter, fixed or estimated
  if(est_phi == 1) {
    for(k in 1:K) {phi_vec[k] = phi[k];}
  } else {
    for(k in 1:K) {phi_vec[k] = 1;}
  }

  // theta is the ma(1) parameter, fixed or estimated
  if(est_theta == 1) {
    for(k in 1:K) {theta_vec[k] = theta[k];}
  } else {
    for(k in 1:K) {theta_vec[k] = 0;}
  }

  for(p in 1:P) {
    sigma_vec[p] = sigma[varIndx[p]]; // convert estimated sigmas to vec form
  }

  if(long_format==0) {
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
  }

  if(proportional_model == 0) {
    for(i in 1:nZ) {
        Z[row_indx[i],col_indx[i]] = z[i]; // convert z to from vec to matrix
      }
      // fill in zero elements in upper diagonal
      if(nZero > 2) {
        for(i in 1:(nZero-2)) {
          Z[row_indx_z[i],col_indx_z[i]] = 0;
        }
      }
      for(k in 1:K) {
        Z[k,k] = zpos[k];// add constraint for Z diagonal
      }
      // this block is for the expansion prior
      for(k in 1:K) {
        if(zpos[k] < 0) {
          indicator[k] = -1;
        } else {
          indicator[k] = 1;
        }
        psi_root[k] = sqrt(psi[k]);
        for(p in 1:P) {
          Z[p,k] = Z[p,k] * indicator[k] * (1/psi_root[k]);
        }
      }
      // initial state for each trend
      for(k in 1:K) {
        x[k,1] = x0[k];
        // trend is modeled as random walk, with optional
        // AR(1) component = phi, and optional MA(1) component
        // theta. Theta is included in the model block below.
        for(t in 2:N) {
          x[k,t] = phi_vec[k]*x[k,t-1] + devs[k,t-1];
        }
      }
      // this block also for the expansion prior, used to convert trends
      for(k in 1:K) {
        //  x[k,1:N] = x[k,1:N] * indicator[k] * psi_root[k];
        for(t in 1:N) {
          x[k,t] = x[k,t] * indicator[k] * psi_root[k];
        }
      }

   }
  if(proportional_model == 1) {
    // initial state for each trend
    for(k in 1:K) {
      x[k,1] = x0[k];
      // trend is modeled as random walk, with optional
      // AR(1) component = phi, and optional MA(1) component
      // theta. Theta is included in the model block below.
      for(t in 2:N) {
        x[k,t] = phi_vec[k]*x[k,t-1] + devs[k,t-1];
      }
    }
    // then try proportional model
    for(p in 1:P) {
      for(k in 1:K) {
        Z[p,k] = p_z[p,k]; // compositions sum to 1 for a time series
      }
    }
  }

  // adjust predictions if process covariates exist
  if(num_pro_covar > 0) {
    for(i in 1:num_pro_covar) {
      // indexed by time, trend, covariate #, covariate value
      x[pro_covar_index[i,2],pro_covar_index[i,1]] = x[pro_covar_index[i,2],pro_covar_index[i,1]] + b_pro[pro_covar_index[i,3], pro_covar_index[i,2]] * pro_covar_value[i];
    }
  }

  // N is sample size, P = time series, K = number trends
  // [PxN] = [PxK] * [KxN]
  pred = Z * x;

  // adjust predictions if observation covariates exist
  if(num_obs_covar > 0) {
    for(i in 1:num_obs_covar) {
      // indexed by time, trend, covariate #, covariate value
      pred[obs_covar_index[i,2],obs_covar_index[i,1]] = pred[obs_covar_index[i,2],obs_covar_index[i,1]] + b_obs[obs_covar_index[i,3], obs_covar_index[i,2]] * obs_covar_value[i];
    }
  }

  if(long_format==1 && est_cor==1) {
    // this is a pain, but we need to calculate the deviations (basically mean y - E[y] for each time point/time series)
    for(n in 1:N) {
      for(p in 1:P) {
        temp_sums[p,n] = 0.0;
        temp_counts[p,n] = 0.0;
      }
    }
    for(i in 1:n_pos) {
      temp_sums[row_indx_pos[i],col_indx_pos[i]] = temp_sums[row_indx_pos[i],col_indx_pos[i]] + (y[i] - pred[row_indx_pos[i],col_indx_pos[i]]);//PxN
      temp_counts[row_indx_pos[i],col_indx_pos[i]] = temp_counts[row_indx_pos[i],col_indx_pos[i]] + 1;
    }
    for(n in 1:N) {
      for(p in 1:P) {
        // Now temp_sums will hold the mean residual for each time - timeseries combination
        temp_sums[p,n] = temp_sums[p,n]/temp_counts[p,n];
      }
    }
    //Omega_derived = multiply_lower_tri_self_transpose(Lcorr);
    Sigma_derived = quad_form_diag(multiply_lower_tri_self_transpose(Lcorr), sigma_vec);

    for(p in 1:P) {
      sigma11 = Sigma_derived[p,p]; //
      Sigma_temp = inverse(subset(Sigma_derived, p, P)); // this is sigma22^-1
      sigma12_vec = subsetvec(Sigma_derived, p, P); // P-1 x 1 matrix
      // conditional mean for multivariate normal, e.g. https://en.wikipedia.org/wiki/Multivariate_normal_distribution
      cond_mean_vec[p] = to_row_vector(sigma12_vec) * Sigma_temp * to_vector(subsetvec2(col(temp_sums,p), p, P));
      // conditional variance of multivariate normal, e.g. https://en.wikipedia.org/wiki/Multivariate_normal_distribution
      cond_sigma_vec[p] = sqrt(sigma11 - to_row_vector(sigma12_vec) * Sigma_temp * to_vector(sigma12_vec));
    }
  }

}
model {
  vector[K] sigma_pro;
  // initial state for each trend
  x0 ~ normal(0, 1); // initial state estimate at t=1
  psi ~ gamma(2, 1); // expansion parameter for par-expanded priors

  // block for process errors - can be estimated or not, and shared or not
  for(k in 1:K) {
    sigma_pro[k] = 1; // default constraint of all DFAs
    if(est_sigma_process==1) {
      if(n_sigma_process==1) {
        sigma_pro[k] = sigma_process[1];
      } else {
        sigma_pro[k] = sigma_process[k];
      }
    }
  }

  // This is deviations - either normal or Student t, and
  // if Student-t, df parameter nu can be estimated or fixed
  for(k in 1:K) {
    if(use_normal == 0) {
      for(t in 1:1) {
        if (estimate_nu == 1) {
          devs[k,t] ~ student_t(nu[1], 0, sigma_pro[k]); // random walk
        } else {
          devs[k,t] ~ student_t(nu_fixed, 0, sigma_pro[k]); // random walk
        }
      }
      for(t in 2:(N-1)) {
        // if MA is not included, theta_vec = 0
        if (estimate_nu == 1) {
          devs[k,t] ~ student_t(nu[1], theta_vec[k]*devs[k,t-1], sigma_pro[k]); // random walk
        } else {
          devs[k,t] ~ student_t(nu_fixed, theta_vec[k]*devs[k,t-1], sigma_pro[k]); // random walk
        }
      }

    } else {
      devs[k,1] ~ normal(0, 1);
      for(t in 2:(N-1)) {
        // if MA is not included, theta_vec = 0
        devs[k,t] ~ normal(theta_vec[k]*devs[k,t-1], sigma_pro[k]);
      }
    }

  }

  // prior for df parameter for t-distribution
  if (estimate_nu == 1) {
    nu[1] ~ gamma(2, 0.1);
  }
  // prior on AR(1) component if included
  if(est_phi == 1) {
    phi ~ normal(0,1); // K elements
  }
  // prior on MA(1) component if included
  if(est_theta == 1) {
    theta ~ normal(0,1); // K elements
  }
  if(est_sigma_process) {
    sigma_process ~ student_t(3, 0, 5);
  }

  if(proportional_model == 0) {
    // prior on loadings
    z ~ normal(0, 1);
    zpos ~ normal(0, 1);// diagonal
  } else {
    for(p in 1:P) {
      p_z[p] ~ dirichlet(alpha_vec);
    }
  }

  // observation variance
  sigma ~ student_t(3, 0, 5);
  if(est_cor == 1) {
    Lcorr ~ lkj_corr_cholesky(1);
  }

  // likelihood for independent
  if(est_cor == 0) {
    if(long_format==0) {
      for(i in 1:P){
        target += normal_lpdf(yall[i] | pred[i], sigma_vec[i]);
      }
    } else {
      for(i in 1:n_pos) {
        //row_indx_pos[i] is the time series, col_index_pos is the time
        target += normal_lpdf(y[i] | pred[row_indx_pos[i],col_indx_pos[i]], sigma_vec[row_indx_pos[i]]);
      }
    }
  } else {
    // need to loop over time slices / columns - each ~ MVN
    if(long_format==0) {
      for(i in 1:N) {
        target += multi_normal_cholesky_lpdf(col(yall,i) | col(pred,i), diag_pre_multiply(sigma_vec, Lcorr));
      }
    } else {
      for(i in 1:n_pos) {
      //  //row_indx_pos[i] is the time series, col_index_pos is the time. note that cond_sigma_vec and cond_mean_vec
      //  // used here, calculated for univariation conditionals above
        target += normal_lpdf(y[i] | pred[row_indx_pos[i],col_indx_pos[i]] + cond_mean_vec[row_indx_pos[i]], cond_sigma_vec[row_indx_pos[i]]);
      }
    }
  }
}
generated quantities {
  vector[n_loglik] log_lik;
  matrix[n_pcor, n_pcor] Omega;
  matrix[n_pcor, n_pcor] Sigma;
  matrix[K,1] xstar; //random walk-trends in future
  vector[K] sigma_pro;
  vector[K] future_devs; // deviations in future
  int<lower=0> j;
  j = 0;

  if(est_cor == 1) {
    Omega = multiply_lower_tri_self_transpose(Lcorr);
    Sigma = quad_form_diag(Omega, sigma_vec);
  }

  // calculate pointwise log_lik for loo package:
  if(est_cor == 0) {
    if(long_format==0) {
      j = 0;
      for(n in 1:N) {
        for(p in 1:P) {
          j = j + 1;
          log_lik[j] = normal_lpdf(yall[p,n] | pred[p,n], sigma_vec[p]);
        }
      }
    } else {
      for(i in 1:n_pos) {
        //row_indx_pos[i] is the time series, col_index_pos is the time
        log_lik[i] = normal_lpdf(y[i] | pred[row_indx_pos[i],col_indx_pos[i]], sigma_vec[row_indx_pos[i]]);
      }
    }

  } else {
    if(long_format==0) {
      for(i in 1:N) {
        log_lik[i] = multi_normal_cholesky_lpdf(col(yall,i) | col(pred,i), diag_pre_multiply(sigma_vec, Lcorr));
      }
    } else {
      for(i in 1:n_pos) {
      //  //row_indx_pos[i] is the time series, col_index_pos is the time
        log_lik[i] = normal_lpdf(y[i] | pred[row_indx_pos[i],col_indx_pos[i]] + cond_mean_vec[row_indx_pos[i]], cond_sigma_vec[row_indx_pos[i]]);
      }
    }
  }

  // future deviations
  // block for process errors - can be estimated or not, and shared or not
  for(k in 1:K) {
    sigma_pro[k] = 1; // default constraint of all DFAs
    if(est_sigma_process==1) {
      if(n_sigma_process==1) {
        sigma_pro[k] = sigma_process[1];
      } else {
        sigma_pro[k] = sigma_process[k];
      }
    }
  }
  for(k in 1:K) {
    if(use_normal == 0) {
      // if MA is not included, theta_vec = 0
      if (estimate_nu == 1) {
        future_devs[k] = student_t_rng(nu[1], theta_vec[k]*devs[k,N-1], sigma_pro[k]); // random walk
      } else {
        future_devs[k] = student_t_rng(nu_fixed, theta_vec[k]*devs[k,N-1], sigma_pro[k]); // random walk
      }
    } else {
      // if MA is not included, theta_vec = 0
      future_devs[k] = normal_rng(theta_vec[k]*devs[k,N-1], sigma_pro[k]);
    }
    xstar[k,1] = x[k,N] + future_devs[k]; // future value of trend at t+1
  }

}
