#' Apply cross validation to DFA model
#'
#' @param stanfit A stanfit object, to preserve the model structure from a call to fit_dfa()
#' @param cv_method The method used for cross validation. The options are 'loocv', where time is ignored and each data point is
#' assigned randomly to a fold. The method 'ltocv' is leave time out cross validation, and time slices are iteratively held out
#' out. Finally the method 'lfocv' implements leave future out cross validation to do one-step ahead predictions.
#' @param fold_ids A vector whose length is the same as the number of total data points. Elements are the fold id of each data point. If not all data points are
#' used (e.g. the lfocv or ltocv approach might only use 10 time steps) the value can be something other than a numbber,
#' e.g. NA
#' @param n_folds Number of folds, defaults to 10
#' @param iter Number of iterations in Stan sampling, defaults to 2000.
#' @param thin Thinning rate in Stan sampling, defaults to 1.
#' @param chains Number of chains in Stan sampling, defaults to 4.
#' @param ... Any other arguments to pass to [rstan::sampling()].
#'
#' @importFrom stats dnorm var
#' @export
#'
#' @examples
#' set.seed(42)
#' s <- sim_dfa(num_trends = 1, num_years = 20, num_ts = 3)
#' obs <- c(s$y_sim[1,], s$y_sim[2,], s$y_sim[3,])
#' long = data.frame("obs" = obs, "ts" = sort(rep(1:3,20)), "time" = rep(1:20,3))
#' m <- fit_dfa(y = long, iter = 50, chains = 1, data_shape="long", sample=FALSE)
#' # random folds
#' fit_cv = dfa_cv(m, cv_method="loocv", n_folds = 5, iter=50, chains=1)
#'
#' # folds can also be passed in
#' fold_ids = sample(1:5, size=nrow(long), replace=TRUE)
#' m <- fit_dfa(y = long, iter = 50, chains = 1, data_shape="long", sample=FALSE)
#' fit_cv = dfa_cv(m, cv_method="loocv", n_folds = 5, iter=50, chains=1, fold_ids=fold_ids)
#'
#' # do an example of leave-time-out cross validation where years are dropped
#' fold_ids = long$time
#' m <- fit_dfa(y = long, iter = 50, chains = 1, data_shape="long", sample=FALSE)
#' fit_cv = dfa_cv(m, cv_method="loocv", iter=100, chains=1, fold_ids = fold_ids)
#' obs_covar = expand.grid("time"=1:20,"timeseries"=1:3,"covariate"=1:2)
#' obs_covar$value=rnorm(nrow(obs_covar),0,0.1)
#'
#' # add in observed data
#' obs <- c(s$y_sim[1,], s$y_sim[2,], s$y_sim[3,])
#' m <- fit_dfa(y = long, iter = 50, chains = 1, obs_covar=obs_covar,data_shape="long", sample=FALSE)
#' fit_cv = dfa_cv(m, cv_method="loocv", n_folds = 5, iter=50, chains=1)
#'
dfa_cv <- function(stanfit,
  cv_method = c("loocv","lfocv"),
  fold_ids = NULL,
  n_folds = 10,
  iter = 2000,
  chains = 4,
  thin = 1,
  ...) {

  cv_method <- match.arg(cv_method, c("loocv","lfocv"))
  if(is.null(fold_ids)) {
    warning("the vector fold_ids containing fold ids is null, so random folds are being used")
    fold_ids <- sample(1:n_folds, nrow(stanfit$orig_data), replace=TRUE)
  }
  if(length(fold_ids) != nrow(stanfit$orig_data)) {
    stop("The length of the vector fold_ids needs to tbe the same as the number of rows in the long format dataframe")
  }
  if(stanfit$shape!="long") {
    stop("Error, please reshape the data into long format")
  }

  if(!is.null(fold_ids)) n_folds = max(fold_ids)
  y <- stanfit$orig_data
  y$time = y$time - min(y$time) + 1

  # loop over the folds, re-fitting the dfa model each time with the folds held out
  log_lik <- matrix(0, nrow=chains*iter/2, ncol = n_folds)
  for(f in 1:n_folds) {

    # fit model holding out each time slice. subset observed data and covar
    y_train <- y
    y_train[which(fold_ids == f),"obs"] <- NA
    y_test <- y[which(fold_ids == f),]
    obs_covar_train=NULL
    if(length(stanfit$sampling_args$data$obs_covar_value) > 0) {
      stanfit$obs_covar$time_timeseries <- paste(stanfit$obs_covar$time,stanfit$obs_covar$timeseries)
      y_train$time_timeseries <- paste(y_train$time,y_train$ts)
      y_test$time_timeseries <- paste(y_test$time,y_test$ts)
      obs_covar_train <- stanfit$obs_covar[which(stanfit$obs_covar$time_timeseries %in% y_train$time_timeseries),1:4]
      obs_covar_test <- stanfit$obs_covar[which(stanfit$obs_covar$time_timeseries %in% y_test$time_timeseries),1:4]
    }
    pro_covar_train=NULL
    if(length(stanfit$sampling_args$data$pro_covar_value) > 0) {
      pro_covar_train <- stanfit$pro_covar[which(fold_ids != f),]
      pro_covar_test <- stanfit$pro_covar[which(fold_ids == f),]
    }

    # fit the new model
    fit.mod <- fit_dfa(y = y_train,
      num_trends = stanfit$sampling_args$data$K,
      varIndx = stanfit$sampling_args$data$varIndx,
      zscore = stanfit$zscore,
      iter = iter,
      chains = chains,
      thin = thin,
      control = stanfit$sampling_args$control,
      nu_fixed = stanfit$sampling_args$data$nu_fixed,
      est_correlation = stanfit$sampling_args$data$est_cor,
      estimate_nu = stanfit$sampling_args$data$estimate_nu,
      estimate_trend_ar = ifelse(stanfit$sampling_args$data$est_phi==1, TRUE, FALSE),
      estimate_trend_ma = ifelse(stanfit$sampling_args$data$est_theta == 1, TRUE, FALSE),
      estimate_process_sigma = ifelse(stanfit$sampling_args$data$est_sigma_process == 1, TRUE, FALSE),
      equal_process_sigma = ifelse(stanfit$sampling_args$data$n_sigma_process == 1, TRUE, FALSE),
      sample = TRUE,
      data_shape = stanfit$shape,
      obs_covar = obs_covar_train,
      pro_covar = pro_covar_train,
      z_bound = stanfit$z_bound,
      z_model = stanfit$z_model,
      verbose = FALSE)

    # extract posterior parameters for the training set
    pars <- rstan::extract(fit.mod$model)
    r <- rotate_trends(fit.mod)
    # loop over each iterations (mcmc sample)
    for(j in 1:nrow(log_lik)) {

      # determine if covariates are included
      obs_covar_offset = rep(0, nrow(y_test))
      if(is.null(obs_covar_train) & is.null(pro_covar_train)) {
        #pred <- pars$Z[j,,] %*% matrix(pars$x[j,,],nrow=stanfit$sampling_args$data$K)
        pred <- r$Z_rot[j,,] %*% matrix(r$trends[j,,],nrow=stanfit$sampling_args$data$K)
        # subset predictions corresponding to observations
        pred <- pred[cbind(y_test$ts,y_test$time)]
        #pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1)
      }
      if(!is.null(obs_covar_train) & is.null(pro_covar_train)) {
        #pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1) + pars$b_obs[j,,] * obs_covar_test$value
        #pred <- pars$Z[j,,] %*% matrix(pars$x[j,,],nrow=stanfit$sampling_args$data$K)
        pred <- r$Z_rot[j,,] %*% matrix(r$trends[j,,],nrow=stanfit$sampling_args$data$K)
        pred <- pred[cbind(y_test$ts,y_test$time)]
        for(i in 1:max(obs_covar_test$covariate)) {
          indx <- which(obs_covar_test$covariate == i)
          pred <- pred + pars$b_obs[j,i,] * obs_covar_test$value[indx]
        }
      }

      log_lik[j,f] <- sum(dnorm(x = y_test$obs,
        mean = pred,
        sd = pars$sigma[j,stanfit$sampling_args$data$varIndx], log=TRUE), na.rm=T)
      #log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
    }

    # Predictions now vary based on how the cross validation is done, and whether covariates used
    #if(cv_method == "loocv") {
    #}
    #if(cv_method == "lfocv") {
      # for(j in 1:nrow(log_lik)) {
      #   # loop over iterations
      #   if(is.null(obs_covar) & is.null(pro_covar)) {
      #     pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1)
      #   }
      #   if(!is.null(obs_covar) & is.null(pro_covar)) {
      #     pred = pars$Z[j,,] %*% matrix(pars$xstar[j,,],ncol=1) + pars$b_obs[j,,] * covar_test$value
      #   }
      #   log_lik[j,k] = sum(dnorm(x = ytest, mean = pred, sd = pars$sigma[j,varIndx], log=TRUE), na.rm=T)
      # }
    #}

  }

  elpds <- apply(log_lik,2,log_sum_exp)
  elpd <- list("log_lik"=log_lik,
    "elpds" = elpds,
    "elpd_kfold"=sum(elpds),
    "se_elpd_kfold" = sqrt(length(elpds) * var(elpds)))
  return(elpd)
}

# more stable than log(sum(exp(x)))
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

# more stable than log(mean(exp(x)))
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}
