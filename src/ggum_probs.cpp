/* These functions are designed to calculate (ggum)
 * - probabilities
 * - first derviatives of probabilties
 * - second derivatives of probabilities
 */

#include <Rcpp.h>
#include "mupp_utilities.h"
using namespace Rcpp;

double exp_ggum(const double theta,
                const double alpha,
                const double delta,
                double tau,
                int    exp_mult = 1){

  // if mult is greater than 3, set tau to 0
  // if mult is outside the range of [1, 3], fix
  if(exp_mult >= 3){
    exp_mult = 3;
    tau      = 0;
  } else if(exp_mult <= 0){
    exp_mult = 1;
  }

  double exp_out = exp(alpha * (exp_mult * (theta - delta) - tau));

  return exp_out;
}

NumericVector q_ggum(const SEXP & thetas,
                     const SEXP & params) {

  // Arguments:
  //  - thetas: a matrix/vector of thetas for one dimension
  //  - params: a matrix of params for one statement [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across one dimensions

  // dimensions of objects
  int n_persons = Rf_xlength(thetas);

  // protect the objects
  PROTECT(thetas);
  PROTECT(params);

  // pointers to objects
  double *pthetas, *pparams;
  pparams = REAL(params);
  pthetas = REAL(thetas);

  //declaring return
  NumericVector probs = no_init(n_persons);

  // declaring parameters
  double alpha = pparams[0];
  double delta = pparams[1];
  double tau   = pparams[2];
  double theta, exp_03, exp_12;

  // calculating probability for each combination ...
  for(int p = 0; p < n_persons; p++){
    theta    = pthetas[p];
    exp_03   = exp_ggum(theta, alpha, delta, tau, 3);
    exp_12   = exp_ggum(theta, alpha, delta, tau, 1) +
               exp_ggum(theta, alpha, delta, tau, 2);
    probs[p] = (1 + exp_03) / (1 + exp_03 + exp_12);
  }

  UNPROTECT(2);

  return probs;

}

//[[Rcpp::export]]
NumericMatrix q_ggum_all(const NumericMatrix & thetas,
                         const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow();
  int n_dims    = params.nrow();
  int n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons);
  NumericVector param = no_init(n_param);
  NumericMatrix probs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta             = thetas.column(dim);
    param             = params.row(dim);
    probs.column(dim) = q_ggum(theta, param);
  }

  return probs;
}

NumericVector pder1_ggum(const SEXP & thetas,
                         const SEXP & params) {

  // Arguments:
  //  - thetas: a matrix/vector of thetas for one dimension
  //  - params: a matrix of params for one statement [alpha, delta, tau]
  // Value:
  //  - P'(Z = 1 | thetas) across one dimensions

  // dimensions of object
  int n_persons = Rf_xlength(thetas);

  // protect the objects
  PROTECT(thetas);
  PROTECT(params);

  // pointers to objects
  double *pthetas, *pparams;
  pparams = REAL(params);
  pthetas = REAL(thetas);

  // declaring return
  NumericVector dprobs = no_init(n_persons);

  // declaring parameters
  double alpha = pparams[0];
  double delta = pparams[1];
  double tau   = pparams[2];
  double theta, exp_1, exp_2, exp_3;

  // calculating derivatives for each combination ...
  for(int p = 0; p < n_persons; p++){
    theta = pthetas[p];
    exp_1 = exp_ggum(theta, alpha, delta, tau, 1);
    exp_2 = exp_ggum(theta, alpha, delta, tau, 2);
    exp_3 = exp_ggum(theta, alpha, delta, tau, 3);
    dprobs[p] = alpha * (exp_1 * (1 - 2 * exp_3) *
                         exp_2 * (2 - 1 * exp_3)) /
                pow((1 + exp_1 + exp_2 + exp_3), 2);
  }

  UNPROTECT(2);

  // (note: derivative of q is -dp)
  return dprobs;
}

// [[Rcpp::export]]
NumericMatrix pder1_ggum_all(const NumericMatrix & thetas,
                             const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P'(Z = 1 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow();
  int n_dims    = params.nrow();
  int n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons);
  NumericVector param = no_init(n_param);
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta              = thetas.column(dim);
    param              = params.row(dim);
    dprobs.column(dim) = pder1_ggum(theta, param);
  }

  return dprobs;
}
