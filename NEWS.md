# mupp 0.0.2

* Updated LICENSE to MIT to make publishing on Github/CRAN easier in the long
  run.
* Updated infrastructure to: 1) add Github workflows for "check-release", "pkgdown",
  and "test-coverage"; 2) add lifecycle stages for the overall package; 3) add
  code coverage; and 4) include badgets for CRAN check, lifecycle, and code
  coverage.
* Fixed overall documentation and standardized package documentation.
* Fixed test cases by adding seeds to reduce random noise, and fixed some style 
  issues.
* `cpp_profiler_tools.cpp` removed due to compilation issues on remote
  systems.

# mupp 0.0.1

* Created functions that simulate responses to the MUPP model, estimate item
  parameters (without constraining any parameters and always using MCMC),
  and estimating person parameters (using either MCMC or EM with BFGS).
* Added C++ functions (using Rcpp) to generate MUPP/GGUM probabilities,
  likelihoods, and corresponding second derivatives.
* Added a test case suite to check simulations and probabilities and ensure
  that they never change with package updates.
