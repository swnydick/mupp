
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mupp

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/swnydick/mupp/workflows/R-CMD-check/badge.svg)](https://github.com/swnydick/mupp/actions)
[![Codecov test
coverage](https://codecov.io/gh/swnydick/mupp/branch/master/graph/badge.svg)](https://codecov.io/gh/swnydick/mupp?branch=master)
<!-- badges: end -->

### An R Packages for the Multi-Unidimensional Pairwise Preference Model

The mupp package contains functions designed to simulation data that
represent the GGUM-based MUPP model with corresponding MCMC and EM
estimation algorithms. Note that estimation for the MUPP model is
assumed to currently be joint estimation of GGUM parameters in a MUPP
framework (which could be problematic and lead to some measurement
error). Currently, there are no options for fixing initial item
parameter estimates or treating items as pretest items, but these
functionalities might come in the future.

## Installation

Within an R session, type the following to install the most current
(development) source package from github

``` r
devtools::install_github("swnydick/mupp")
```

## Help

Take a look at the help pages for `simulate_mupp_parameters`,
`simulate_mupp_responses`, `estimate_mupp_parameters`, and
`estimate_mupp_thetas` to get started.

The package homepage is currently located
[here](https://swnydick.github.io/mupp/index.html).
