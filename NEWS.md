# bayesdfa 0.1.0

* Initial submission to CRAN.

# bayesdfa 0.1.1

* Changed Makevars per exchange with Stan developers.

# bayesdfa 0.1.2

* Changed find_inverted_chains() and invert_chains() to be compatible with dplyr 0.8 release. Specifically, removed deprecated group_by_() and summarise_() functions and changed code to remove unused factor levels.

# bayesdfa 0.1.3

* Changed Makevars and Makevars.win to rely on CXX_STD = CXX14 from CXX11. Also added vignette related to inclusion of covariates

# bayesdfa 0.1.5

* Added additional functionality to relax limits on AR(1) parameter (phi), MA(1) parameter (theta), and flexibility in estimated the standard deviation of latent trends. Also modified the data object passed in to be either a wide matrix (as previously done) or a long format data frame. The latter allows for multiple observations / time step. Finally, an additional and alternative constraint was introduced for Z, allowing elements to be modeled as a Dirchlet process, rather than conventional DFA. 

# bayesdfa 0.1.6

* Removed warning related to vignette and noLD test

# bayesdfa 0.1.7

* Added non-gaussian families (poisson, negative binomial, bernoulli, Gamma, lognormal). Also included a function for doing cross validation and calculating the expected log posterior density. Another new feature included smooth models (Gaussian process, B-splines) as alternative models for trends conventionally modeled as random walks. Added functions dfa_trends(), dfa_loadings() and dfa_fitted() for extracting trends, loadings, and fitted values. 

# bayesdfa 1.0.0

* Added constraint on diagonal of Z matrix to keep parameter estimates from 'flipping' within MCMC chains. Ensures convergence for problematic cases. This was present in 0.1.1, but later removed. 

# bayesdfa 1.1.0

* Following 1.0.0, included a new argument to fit_dfa() function 'expansion_prior' that allows user to toggle on / off the constraint. If not included (default=FALSE), there is no constraint on the Z diagonal, and post-hoc MCMC chain inverting resolves identifiability. If 'expansion_prior' = TRUE, then the positive constraint is applied, in combination with the expansion prior for trends and loadings. 

# bayesdfa 1.2.0

Add penalized spline models, so that the 'trend_model' argument may take on
"rw" for conventional random walks, "bs" for B-splines, "ps" for "P-splines",
or "gp" for Gaussian processes

# bayesdfa 1.3.0

Change to new Stan syntax

# bayesdfa 1.3.1

Versioning

# bayesdfa 1.3.2

- Add compatibility with new rstan 
- Changed weights argument to 'inv_var_weights' and 'likelihood_weights' for the glmmTMB/sdmTMB/brms style
