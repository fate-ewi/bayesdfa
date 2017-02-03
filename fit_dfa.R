fit_dfa = function(y = y,
  num_trends = 2,
  varIndx = NULL,
  zscore = TRUE,
  iter = 4000,
  chains = 1,
  control = list(adapt_delta = 0.99),
  nu = 7,
  model = c("dfa.stan", "tvdfa.stan")) {
  # parameters for DFA
  N = ncol(y)
  P = nrow(y)
  K = num_trends # number of dfa trends
  nZ = P * K - sum(1:K) + K
  
  if (zscore == TRUE) {
    for (i in 1:P) {
      y[i, ] = scale(y[i, ], center = TRUE, scale = TRUE)
    }
  }
  
  mat_indx = matrix(0, P, K)
  start = 1
  for (k in 1:K) {
    if (k == 1)
      mat_indx[, k] = (1:nZ)[start:(start + P - k)]
    if (k > 1)
      mat_indx[-c(0:(k - 1)), k] = (1:nZ)[start:(start + P - k)]
    start = start + (P - k + 1)
  }
  
  row_indx = matrix((rep(1:P, K)), P, K)[which(mat_indx > 0)]
  col_indx = rep(1:K, times = P:(P - K + 1))
  row_indx_z = matrix((rep(1:P, K)), P, K)[which(mat_indx == 0)]
  col_indx_z = matrix(sort(rep(1:K, P)), P, K)[which(mat_indx == 0)]
  row_indx_z = c(row_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  col_indx_z = c(col_indx_z, 0, 0)# +2 zeros for making stan ok with data types
  nZero = length(row_indx_z)
  
  # set the model up to have shared variances between first two time series,
  # third is different
  if (is.null(varIndx))
    varIndx = rep(1, P)
  nVariances = length(unique(varIndx))
  
  # indices of positive values - stan can't handle NAs
  row_indx_pos = matrix((rep(1:P, N)), P, N)[which(!is.na(y))]
  col_indx_pos = matrix(sort(rep(1:N, P)), P, N)[which(!is.na(y))]
  n_pos = length(row_indx_pos)
  y = y[which(!is.na(y))]
  
  data_list = list(
    N,
    P,
    K,
    nZ,
    y,
    row_indx,
    col_indx,
    nZero,
    row_indx_z,
    col_indx_z,
    nZero,
    row_indx_z,
    col_indx_z,
    row_indx_pos,
    col_indx_pos,
    n_pos,
    nu
  )
  pars <- c("x", "Z", "sigma", "log_lik")
  if (model[[1]] == "tvdfa.stan") pars <- c(pars, "tau")
  mod = stan(
    data = data_list,
    pars = pars,
    file = model[[1]],
    chains = chains,
    iter = iter,
    control = control
  )
}
