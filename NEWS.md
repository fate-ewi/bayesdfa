# bayesdfa 0.1.0

* Initial submission to CRAN.

# bayesdfa 0.1.1

* Changed Makevars per exchange with Stan developers.

# bayesdfa 0.1.2

* Changed find_inverted_chains() and invert_chains() to be compatible with dplyr 0.8 release. Specifically, removed deprecated group_by_() and summarise_() functions and changed code to remove unused factor levels.

# bayesdfa 0.1.3

* Changed Makevars and Makevars.win to rely on CXX_STD = CXX14 from CXX11. Also added vignette related to inclusion of covariates
