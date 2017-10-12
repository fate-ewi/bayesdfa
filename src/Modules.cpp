#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4corr_mod) {


    class_<rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> >("model_corr")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_corr_namespace::model_corr, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4dfa_mod) {


    class_<rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> >("model_dfa")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_dfa_namespace::model_dfa, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4regime_1_mod) {


    class_<rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> >("model_regime_1")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_regime_1_namespace::model_regime_1, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4regime_2_mod) {


    class_<rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> >("model_regime_2")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_regime_2_namespace::model_regime_2, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4regime_3plus_mod) {


    class_<rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> >("model_regime_3plus")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_regime_3plus_namespace::model_regime_3plus, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
#include <Rcpp.h>
using namespace Rcpp ;
#include "include/models.hpp"

RCPP_MODULE(stan_fit4tvdfa_fixed_mod) {


    class_<rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> >("model_tvdfa_fixed")

    .constructor<SEXP,SEXP>()


    .method("call_sampler", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::call_sampler)
    .method("param_names", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_names)
    .method("param_names_oi", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_names_oi)
    .method("param_fnames_oi", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_fnames_oi)
    .method("param_dims", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_dims)
    .method("param_dims_oi", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_dims_oi)
    .method("update_param_oi", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::update_param_oi)
    .method("param_oi_tidx", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::param_oi_tidx)
    .method("grad_log_prob", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::grad_log_prob)
    .method("log_prob", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::log_prob)
    .method("unconstrain_pars", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::unconstrain_pars)
    .method("constrain_pars", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::constrain_pars)
    .method("num_pars_unconstrained", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::num_pars_unconstrained)
    .method("unconstrained_param_names", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::unconstrained_param_names)
    .method("constrained_param_names", &rstan::stan_fit<model_tvdfa_fixed_namespace::model_tvdfa_fixed, boost::random::ecuyer1988> ::constrained_param_names)
    ;
}
