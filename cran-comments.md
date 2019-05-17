## Test environments

* local OS X install, R 3.5.1
* ubuntu 14.04.5 LTSR (on travis-ci), R Under development (unstable)
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 2 notes

* installed size is  5.7Mb; sub-directories of 1Mb or more: libs 4.9Mb
* GNU make is a SystemRequirements

Explanation: this is from the compiled 'Stan' model and associated libraries, and is necessary https://mc-stan.org/rstantools/articles/minimal-rstan-package.html.

Note that we have an associated manuscript accepted in The R Journal's next issue.
