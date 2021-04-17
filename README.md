# mupp

### An R Packages for the Multi-Unidimensional Pairwise Preference Model

The mupp package contains functions designed to simulation data that represent
the GGUM-based MUPP model with corresponding MCMC and EM estimation algorithms.
Note that estimation for the MUPP model is assumed to currently be joint
estimation of GGUM parameters in a MUPP framework (which could be problematic
and lead to some measurement error). Currently, there are no options for
fixing initial item parameter estimates or treating items as pretest items,
but these functionalities might come in the future.

Within an R session, type the following to install the most current (development)
source package from github

``` r
if(!require(devtools)){
  install.packages("devtools")
}
devtools::install_github("swnydick/mupp")
```

Take a look at the help pages for `simulate_mupp_parameters`,
`simulate_mupp_responses`, `estimate_mupp_parameters`, and
`estimate_mupp_thetas` to get started.
