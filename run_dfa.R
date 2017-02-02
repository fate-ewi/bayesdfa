library(rstan)

y = matrix(runif(10*3), 3, 10) # 3 time series, 10 years

fit_dfa = function(y = y, num_trends = 2, varIndx = NULL) {
  
  # parameters for DFA
  N = ncol(y)
  P = nrow(y)
  K = num_trends # number of dfa trends
  nZ = P*K - sum(1:K) + K
  
  mat_indx = matrix(0, P, K)
  start = 1
  for(k in 1:K) {
    if(k==1) mat_indx[,k] = (1:nZ)[start:(start+P-k)] 
    if(k > 1) mat_indx[-c(0:(k-1)),k] = (1:nZ)[start:(start+P-k)] 
    start = start + (P - k + 1)
  }
  
  row_indx = matrix((rep(1:P, K)), P, K)[which(mat_indx > 0)]
  col_indx = rep(1:K, times = P:(P-K+1))
  row_indx_z = matrix((rep(1:P, K)), P, K)[which(mat_indx == 0)]
  col_indx_z = matrix(sort(rep(1:K, P)), P, K)[which(mat_indx == 0)]
  nZero = length(row_indx_z)
  
  # set the model up to have shared variances between first two time series,
  # third is different
  if(is.null(varIndx)) varIndx = rep(1,P)
  nVariances = length(unique(varIndx))
  
  data_list = list("N","P","K","nZ","y","row_indx",
    "col_indx","nZero","row_indx_z","col_indx_z","nZero",
    "row_indx_z", "col_indx_z")
  mod = stan(data = data_list, 
    pars = c("Z","sigma"), file="dfa.stan", 
    chains = 1, iter=500, thin=1)
}



