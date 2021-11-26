## Test environments

* local OS X install, R 3.6.0 (and more recent)
* ubuntu 14.04.5 LTSR (on travis-ci), R Under development (unstable)
* win-builder (devel and release)

## R CMD check results

0 errors | 1 warnings | 2 notes

* installed size is  7.3Mb; sub-directories of 1Mb or more: libs 5.8Mb
* GNU make is a SystemRequirements
* Warning A complete check needs the 'checkbashisms' script

Explanation: this is from the compiled 'Stan' model and associated libraries, and is necessary https://mc-stan.org/rstantools/articles/minimal-rstan-package.html.

Note that we have an associated manuscript accepted in The R Journal. Citation: 
Eric J. Ward, Sean C. Anderson, Luis A. Damiano, Mary E. Hunsicker and Michael A. Litzow , The R Journal (2019) 11:2, pages 46-55, https://journal.r-project.org/archive/2019/RJ-2019-007/RJ-2019-007.pdf
