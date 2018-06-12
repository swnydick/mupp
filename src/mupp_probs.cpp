/* These functions are designed to calculate (mupp)
 * - probabilities
 * - first derviatives of probabilties
 * - second derivatives of probabilities
 */

/* _impl stands for either:
 * - implicit (these functions are implicitly called in R)
 * - implementation (these functions are the implementation in C++)
 * - this follows Hadley convention so that no R file overwrites a C++ impl.
 */

#include <Rcpp.h>
using namespace Rcpp;

// /* HELPER FUNCTIONS
//  *   - select_cols (select non-proximate columns)
//  */  - select_rows (select non-prosimate rows)
NumericMatrix select_cols(NumericMatrix X,
                          IntegerVector ind){

 // Arguments:
 //  - X:   a matrix to select columns
 //  - ind: an integer vector of columns to select

 // capturing constants
 int out_cols = ind.size();

 // declaring output
 NumericMatrix Y(X.nrow(), out_cols);

 // iterating
 for(int i = 0; i < out_cols; i++){
   Y(_, i) = X(_, ind[i]);
 }

  return Y;
}

NumericMatrix select_rows(NumericMatrix X,
                          IntegerVector ind){

 // Arguments:
 //  - X:   a matrix to select columns
 //  - ind: an integer vector of columns to select

 // capturing constants
 int out_rows = ind.size();

 // declaring output
 NumericMatrix Y(out_rows, X.cols());

 // iterating
 for(int i = 0; i < out_rows; i++){
   Y(i, _) = X(ind[i], _);
 }

  return Y;
}


/* GGUM STUFF
 *  - exp_ggum
 *  - p_ggum/q_ggum_all
 *  - p_der1_theta_ggum
 */
NumericVector exp_ggum(NumericVector thetas,
                       NumericVector params,
                       int exp_mult = 1){

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - exp(alpha(theta - delta) - tau)

  // declaring parameters
  double alpha = params[0];
  double delta = params[1];
  double tau   = params[2];

  // if mult is greater than 3, set tau to 0
  // if mult is outside the range of [1, 3], fix
  if(exp_mult >= 3){
    exp_mult = 3;
    tau      = 0;
  } else if(exp_mult <= 0){
    exp_mult = 1;
  }

  NumericVector exp_1 = exp(alpha * (exp_mult * (thetas - delta) - tau));

  return exp_1;
}

// [[Rcpp::export]]
NumericVector p_ggum(NumericVector thetas,
                     NumericVector params){

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z_s = 1 | thetas)

  // calculate exponent expressions
  NumericVector exp_12 = exp_ggum(thetas, params, 1) +
                         exp_ggum(thetas, params, 2);
  NumericVector exp_3  = exp_ggum(thetas, params, 3);

  // calculate probabilities
  NumericVector probs  = exp_12 / (1 + exp_12 + exp_3);

  return probs;
}

// [[Rcpp::export]]
NumericVector q_ggum_all(NumericMatrix thetas,
                         NumericMatrix params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow();
  int n_dims    = params.nrow();

  // indicating return vector
  NumericMatrix probs(n_persons, n_dims);

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    probs(_, dim) = 1 - p_ggum(thetas(_, dim), params(dim, _));
  }

  return probs;
}

NumericVector pder1_theta_ggum(NumericVector thetas,
                               NumericVector params){

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P'(Z_s = 1 | thetas)

  // declaring parameters
  double alpha = params[0];

  // calculate exponent expressions
  NumericVector exp_1 = exp_ggum(thetas, params, 1);
  NumericVector exp_2 = exp_ggum(thetas, params, 2);
  NumericVector exp_3 = exp_ggum(thetas, params, 3);

  // calculate derivative of p1 (derivative of p0 is -dp1)
  NumericVector dprobs = alpha * (exp_1 * (1 - 2 * exp_3) +
                                  exp_2 * (2 - 1 * exp_3)) /
                         pow((1 + exp_1 + exp_2 + exp_3), 2);

  return dprobs;
}

/* MUPP STUFF
 *   - p_mupp_pick0 - numerator (individual component) of PICK probability
 *   - p_mupp_pick1 - OVERALL PICK probability
 *   - p_mupp_rank1 - OVERALL RANK probability
 *   - p_mupp_rank_impl - OVERALL RANK probability for an individual order OR
 *                        all possible order combinations
 */

// [[Rcpp::export]]
NumericVector p_mupp_pick0(NumericMatrix Q,
                           int picked_dim = 1){

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - picked_dim: the picked dimension
  // Value:
  //  - P(Z_s = 1 | thetas) individual element where s is the selected option

  // declare number of persons/dimensions
  int n_persons = Q.nrow();
  int n_dims    = Q.ncol();

  // vectors to store stuff
  NumericVector probs(n_persons, 1.0);

  // calculating (intersection) probability across picked dimension
  for(int dim = 0; dim < n_dims; dim++){
    if(dim == picked_dim){
      probs = probs * (1 - Q(_, dim));
    } else{
      probs = probs * Q(_, dim);
    }
  }

  return probs;
}

// [[Rcpp::export]]
NumericVector p_mupp_pick1(NumericMatrix Q,
                           int picked_dim = 1){

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - picked_dim: the picked dimension
  // Value:
  //  - P(Z_s = 1 | thetas) OVERALL where s is the selected option

  // declare number of persons/dimensions
  int n_persons = Q.nrow();
  int n_dims    = Q.ncol();

  // vectors to store stuff
  NumericVector numer(n_persons);
  NumericVector denom(n_persons);
  NumericVector p(n_persons);

  // calculating probability across picked dimension
  for(int dim = 0; dim < n_dims; dim++){

    // determine individual probability multiplication
    p = p_mupp_pick0(Q, dim);

    // add to numerator IF picked_dim, otherwise add to denominator
    if(dim == picked_dim){
      numer = numer + p;
    } else{
      denom = denom + p;
    }
  }

  return numer / (numer + denom);
}

// [[Rcpp::export]]
NumericVector p_mupp_rank1(NumericMatrix Q,
                           IntegerVector order){

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - order: the order of the output
  // Value:
  //  - P(Z_s = 1 | thetas) where s is the selected option

  // declare number of persons/dimensions
  int n_persons = Q.nrow();
  int n_dims    = order.size();

  // matrices to store stuff
  NumericMatrix Q_order = select_cols(Q, order);

  // indicating return vector
  NumericVector probs(n_persons, 1.0);

  // calculating probability of being in particular order (1 by 1)
  for(int dim = 0; dim < n_dims - 1; dim++){
    probs = probs * p_mupp_pick1(Q_order(_, Range(dim, n_dims - 1)), 0);
  }

  return probs;
}

// // [[Rcpp::export]]
// NumericMatrix p_mupp_rank_impl(NumericMatrix thetas,
//                                NumericMatrix params,
//                                IntegerVector dims) {
//   // Arguments:
//   //  - thetas: matrix of persons x dims (for all dims)
//   //  - params: matrix of dims x params (for single item dims)
//   //  - dims:   vector of dims of the items (in sorted order)
//   // Value:
//   //  - matrix of P(s > t > ...)(theta_s, theta_t, theta_...) for all s, t, ... in dims
//
//   // declare number of items and persons
//   int n_items   = std::max(params_s.nrow(), params_t.nrow());
//   int n_persons = std::max(thetas_s.size(), thetas_t.size());
//   int n_params  = std::max(params_s.ncol(), params_t.ncol());
//
//   // indicate temporary storage vectors
//   NumericVector params_si(n_params),
//                 params_ti(n_params),
//                 p_s1(n_persons),
//                 p_t0(n_persons);
//
//   // indicate return matrix
//   NumericMatrix probs(n_persons, n_items);
//
//   // cycle through and add vectors to return matrix
//   for(int item = 0; item < n_items; item++){
//
//     // params for this specific item
//     params_si = params_s(item, _);
//     params_ti = params_t(item, _);
//
//     // probability of s/t given ggum model (can use other models?)
//     p_s1 = p_ggum_impl(thetas_s, params_si);
//     p_t0 = 1 - p_ggum_impl(thetas_t, params_ti);
//
//     // set probs to item column of probability matrix
//     probs(_, item) = (p_s1 * p_t0) / (p_s1 * p_t0 + (1 - p_s1) * (1 - p_t0));
//   }
//
//   return probs;
//
// }

