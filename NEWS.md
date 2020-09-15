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
