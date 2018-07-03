/* These functions are designed to calculate (ggum)
 * - probabilities
 * - first derviatives of probabilties
 * - second derivatives of probabilities
 */

#include <Rcpp.h>
#include "mupp_utilities.h"
using namespace Rcpp;

// HELPER //

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

// PROBABILITY //

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
  double alpha = pparams[0],
         delta = pparams[1],
         tau   = pparams[2];
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

NumericMatrix q_ggum_all(const NumericMatrix & thetas,
                         const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow(),
      n_dims    = params.nrow(),
      n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons),
                param = no_init(n_param);
  NumericMatrix probs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta             = thetas.column(dim);
    param             = params.row(dim);
    probs.column(dim) = q_ggum(theta, param);
  }

  return probs;
}

NumericMatrix p_ggum_all(const NumericMatrix & thetas,
                         const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow(),
      n_dims    = params.nrow(),
      n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons),
                param = no_init(n_param);
  NumericMatrix probs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta             = thetas.column(dim);
    param             = params.row(dim);
    probs.column(dim) = 1 - q_ggum(theta, param);
  }

  return probs;
}

// DERIVATIVE (w.r.t. THETA) //

// numerator of pder1_ggum derivative
double pder1_numer_ggum(double exp_1,
                        double exp_2,
                        double exp_3,
                        double alpha){

  double numer = alpha * (exp_1 * (1 - 2 * exp_3) +
                          exp_2 * (2 - 1 * exp_3));

  return numer;
}

// denominator of pder1_ggum derivative
double pder1_denom_ggum(double exp_1,
                        double exp_2,
                        double exp_3){

  double denom = pow(1 + exp_1 + exp_2 + exp_3, 2);

  return denom;
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
  double alpha = pparams[0],
         delta = pparams[1],
         tau   = pparams[2];
  double theta, exp_1, exp_2, exp_3;

  // calculating derivatives for each combination ...
  for(int p = 0; p < n_persons; p++){
    theta = pthetas[p];
    exp_1     = exp_ggum(theta, alpha, delta, tau, 1);
    exp_2     = exp_ggum(theta, alpha, delta, tau, 2);
    exp_3     = exp_ggum(theta, alpha, delta, tau, 3);
    dprobs[p] = pder1_numer_ggum(exp_1, exp_2, exp_3, alpha) /
                pder1_denom_ggum(exp_1, exp_2, exp_3);
  }

  UNPROTECT(2);

  // (note: derivative of q is -dp)
  return dprobs;
}

NumericMatrix pder1_ggum_all(const NumericMatrix & thetas,
                             const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P'(Z = 1 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow(),
      n_dims    = params.nrow(),
      n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons),
                param = no_init(n_param);
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta              = thetas.column(dim);
    param              = params.row(dim);
    dprobs.column(dim) = pder1_ggum(theta, param);
  }

  return dprobs;
}

// SECOND DERIVATIVE (w.r.t. THETA) //

NumericVector pder2_ggum(const SEXP & thetas,
                         const SEXP & params) {

  // Arguments:
  //  - thetas: a matrix/vector of thetas for one dimension
  //  - params: a matrix of params for one statement [alpha, delta, tau]
  // Value:
  //  - P''(Z = 1 | thetas) across one dimensions

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
  double alpha = pparams[0],
         delta = pparams[1],
         tau   = pparams[2];
  double theta, exp_1, exp_2, exp_3, numer, dnumer, denom, ddenom;

  // calculating second derivatives for each combination ...
  // note: tried to simplify, but the simplification doesn't seem much simpler ...
  for(int p = 0; p < n_persons; p++){
    theta  = pthetas[p];
    exp_1  = exp_ggum(theta, alpha, delta, tau, 1);
    exp_2  = exp_ggum(theta, alpha, delta, tau, 2);
    exp_3  = exp_ggum(theta, alpha, delta, tau, 3);

    // first derivative numerator and denominator
    numer  = pder1_numer_ggum(exp_1, exp_2, exp_3, alpha);
    denom  = pder1_denom_ggum(exp_1, exp_2, exp_3);

    // derivative of first derivative numerator and denominator
    dnumer = pow(alpha, 2) * (exp_1 * (1 - 8 * exp_3) +
                              exp_2 * (4 - 5 * exp_3));
    ddenom = 2 * alpha * (1 + exp_1 + exp_2 + exp_3) * (exp_1 + 2 * exp_2 + 3 * exp_3);

    // standard division rule for a scalar
    dprobs[p] = (denom * dnumer - numer * ddenom) / pow(denom, 2);
  }

  UNPROTECT(2);

  // (note: second derivative of q is -d2p)
  return dprobs;
}

NumericMatrix pder2_ggum_all(const NumericMatrix & thetas,
                             const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P''(Z = 1 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow(),
      n_dims    = params.nrow(),
      n_param   = params.ncol();

  // indicating return vector
  NumericVector theta = no_init(n_persons),
                param = no_init(n_param);
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    theta              = thetas.column(dim);
    param              = params.row(dim);
    dprobs.column(dim) = pder2_ggum(theta, param);
  }

  return dprobs;
}
