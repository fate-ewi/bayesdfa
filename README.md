<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesdfa

[![R build
status](https://github.com/fate-ewi/bayesdfa/workflows/R-CMD-check/badge.svg)](https://github.com/fate-ewi/bayesdfa/actions)

bayesdfa implements Bayesian Dynamic Factor Analysis (DFA) with Stan.

You can install the development version of the package with:

``` r
# install.packages("devtools")
devtools::install_github("fate-ewi/bayesdfa")
```

## Overview

A brief video overview of the package is here,

<figure class="video_container">
<iframe width="560" height="315" src="https://www.youtube.com/embed/yTX7D8_Ad8g" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen>
</iframe>
</figure>

## Vignettes

We’ve put together several vignettes for using the `bayesdfa` package.  
[Overview](https://fate-ewi.github.io/bayesdfa/articles/a1_bayesdfa.html)  
[Combining
data](https://fate-ewi.github.io/bayesdfa/articles/a2_combining_data.html)  
[Including
covariates](https://fate-ewi.github.io/bayesdfa/articles/a3_covariates.html)  
[Smooth trend
models](https://fate-ewi.github.io/bayesdfa/articles/a4_smooth.html)  
[Estimating process
variance](https://fate-ewi.github.io/bayesdfa/articles/a5_estimate_process_sigma.html)  
[Compositional
models](https://fate-ewi.github.io/bayesdfa/articles/a6_compositional.html)  
[DFA for big
data](https://fate-ewi.github.io/bayesdfa/articles/a7_bigdata.html).

Additional examples can be found in the course that Eli Holmes, Mark
Scheuerell, and Eric Ward teach at the University of Washington:  
[Course webpage](https://nwfsc-timeseries.github.io/atsa/)  
[Lab book](https://nwfsc-timeseries.github.io/atsa/)

## Citing

For DFA models in general, we recommend citing the MARSS package or user
guide.

    @article{marss_package,
        title = {{MARSS}: multivariate autoregressive state-space models for analyzing time-series data},
        volume = {4},
        url = {https://pdfs.semanticscholar.org/5d41/b86dff5f977a0eac426a924cf7917220fc9a.pdf},
        number = {1},
        journal = {R Journal},
        author = {Holmes, E.E. and Ward, Eric J. and Wills, K.},
        year = {2012},
        pages = {11--19}
    }

    @article{marss_user_guide,
        title = {{MARSS}: Analysis of multivariate timeseries using the MARSS package},
        url = {https://cran.r-project.org/web/packages/MARSS/vignettes/UserGuide.pdf},
        author = {Holmes, E.E. and Scheurell, M.D. and Ward, Eric J.},
        year = {2020},
    }

For citing the `bayesdfa` package using Bayesian estimation, or models
with extra features (such as extremes), cite

<https://journal.r-project.org/archive/2019/RJ-2019-007/index.html>

    @article{ward_etal_2019,
      author = {Eric J. Ward and Sean C. Anderson and Luis A. Damiano and
              Mary E. Hunsicker and Michael A. Litzow},
      title = {{Modeling regimes with extremes: the bayesdfa package for
              identifying and forecasting common trends and anomalies in
              multivariate time-series data}},
      year = {2019},
      journal = {{The R Journal}},
      doi = {10.32614/RJ-2019-007},
      url = {https://journal.r-project.org/archive/2019/RJ-2019-007/index.html}
    }

### Applications

The ‘bayesdfa’ models were presented to the PFMC’s SSC in November 2017
and have been included in the 2018 California Current Integrated
Ecosystem Report,
<https://www.integratedecosystemassessment.noaa.gov/Assets/iea/california/Report/pdf/CCIEA-status-report-2018.pdf>

### Funding

The ‘bayesdfa’ package was funded by a NOAA Fisheries and the
Environment (FATE) grant on early warning indicators, led by Mary
Hunsicker and Mike Litzow.

### NOAA Disclaimer

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

<img src="https://raw.githubusercontent.com/nmfs-general-modeling-tools/nmfspalette/main/man/figures/noaa-fisheries-rgb-2line-horizontal-small.png" height="75" alt="NOAA Fisheries">

[U.S. Department of Commerce](https://www.commerce.gov/) \| [National
Oceanographic and Atmospheric Administration](https://www.noaa.gov) \|
[NOAA Fisheries](https://www.fisheries.noaa.gov/)
