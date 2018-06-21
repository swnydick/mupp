/* BENCHMARKING FUNCTIONS
 * - These are simply benchmarking functions to try to make Rcpp code faster.
 * - See mupp_probs.cpp for explanations of these functions.
 */

#include <Rcpp.h>
using namespace Rcpp;

///////////////////////////
// FIND ALL PERMUTATIONS //
//////////////////////////

// [[Rcpp::export]]
IntegerMatrix find_all_permutations(int n,
                                    int init = 0){

  // Arguments:
  //  - n: an integer greater than or equal to 1
  // Value:
  //  - permutation matrix in reasonable order

  // Argument Checks
  if(n < 1){
    stop("n must be a scalar greater than 1");
  }

  // declaring helper vectors
  IntegerVector n_r(1, n);

  // declaring parameter parts
  int n_rows   = factorial(n_r)[0];
  int n_cols   = n;
  int i        = 0;

  // indicating return matrix and temporary storage vector
  IntegerVector cur_perm = seq_len(n_cols) - (1 - init);
  IntegerMatrix all_perms(n_rows, n_cols);

  // add permutation to all_perms vector, then
  // use iterators to find next permutation
  do{
    all_perms(i, _) = cur_perm;
    i              += 1;
  } while (std::next_permutation(cur_perm.begin(), cur_perm.end()));

  return all_perms;
}


/////////////////////////////
// SELECT COLS/SELECT ROWS //
////////////////////////////

// [[Rcpp::export]]
NumericMatrix select_cols0(NumericMatrix X,
                           IntegerVector ind) {

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


// [[Rcpp::export]]
NumericMatrix select_cols(NumericMatrix X,
                          IntegerVector ind) {

  // Arguments:
  //  - X:   a matrix to select columns
  //  - ind: an integer vector of columns to select

  // capturing constants
  int out_cols = ind.size();
  int out_rows = X.nrow();

  NumericMatrix Y(Rf_allocMatrix(REALSXP, out_rows, out_cols));

   for(int j = 0; j < out_cols; j++){
     for(int i = 0; i < out_rows; i++){
       Y[i + j * out_rows] = X[i + ind[j] * out_rows];
     }
   }

  return Y;
}

/*** R
library(microbenchmark)
nrow <- 20
ncol <- 4000

mat  <- matrix(rnorm(nrow * ncol), nrow = nrow)
ind  <- sample(x = ncol, size = 10) - 1

microbenchmark(x <- select_cols0(mat, ind),
               y <- select_cols(mat, ind))

*/

NumericMatrix select_rows(NumericMatrix X,
                          IntegerVector ind) {

 // Arguments:
 //  - X:   a matrix to select columns
 //  - ind: an integer vector of columns to select

 // capturing constants
 int out_rows = ind.size();

 // declaring output
 NumericMatrix Y(out_rows, X.ncol());

 // iterating
 for(int i = 0; i < out_rows; i++){
   Y(i, _) = X(ind[i], _);
 }

  return Y;
}

//////////////
// EXP GGUM //
//////////////

// [[Rcpp::export]]
NumericVector exp_ggum0(NumericVector thetas,
                        NumericVector params,
                        int exp_mult = 1) {

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
NumericVector exp_ggum(SEXP thetas,
                       SEXP params,
                       int exp_mult = 1) {

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
  NumericVector exp_1 = no_init(n_persons);

  // declaring parameters
  double alpha = pparams[0];
  double delta = pparams[1];
  double tau   = pparams[2];

  // if mult is greater than 3, set tau to 0
  // if mult is outside the range of [1, 3], fix
  if(exp_mult >= 3){
    exp_mult = 3;
    tau      = 0;
  } else if(exp_mult <= 0){
    exp_mult = 1;
  }

  for(int p = 0; p < n_persons; p++){
    exp_1[p] = exp(alpha * (exp_mult * (pthetas[p] - delta) - tau));
  }

  UNPROTECT(2);

  return exp_1;

}

/*** R
library(microbenchmark)

thetas <- rnorm(100000)
params <- c(1, 0, -1)

microbenchmark(x <- exp_ggum0(thetas, params),
               y <- exp_ggum(thetas, params))

# exp_ggum is faster
*/

////////////
// Q GGUM //
////////////

// [[Rcpp::export]]
NumericVector q_ggum0(NumericVector thetas,
                      NumericVector params) {

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z_s = 0 | thetas) for a single parameter set

  // calculate exponent expressions
  NumericVector exp_12 = exp_ggum(thetas, params, 1) +
                         exp_ggum(thetas, params, 2);
  NumericVector exp_03 = 1 +
                         exp_ggum(thetas, params, 3);

  // calculate probabilities
  NumericVector probs  = exp_03 / (exp_12 + exp_03);

  return probs;
}


// [[Rcpp::export]]
NumericVector q_ggum1(SEXP thetas,
                      SEXP params) {

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z_s = 0 | thetas) for a single parameter set

  // calculate exponent expressions
  NumericVector exp_03 = 1 +
                         exp_ggum(thetas, params, 3);

  // calculate probabilities
  NumericVector probs  = (exp_03) /
                         (exp_03 +
                          exp_ggum(thetas, params, 1) + exp_ggum(thetas, params, 2));

  return probs;
}


// [[Rcpp::export]]
NumericVector q_ggum2(SEXP thetas,
                      SEXP params) {

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

  // // if mult is greater than 3, set tau to 0
  // // if mult is outside the range of [1, 3], fix
  // if(exp_mult >= 3){
  //   exp_mult = 3;
  //   tau      = 0;
  // } else if(exp_mult <= 0){
  //   exp_mult = 1;
  // }

  double exp_03, exp_12;
  for(int p = 0; p < n_persons; p++){
    exp_03   = exp(alpha * (3 * (pthetas[p] - delta)));
    exp_12   = exp(alpha * (1 * (pthetas[p] - delta) - tau)) +
               exp(alpha * (2 * (pthetas[p] - delta) - tau));
    probs[p] = (1 + exp_03) / (1 + exp_03 + exp_12);
  }

  UNPROTECT(2);

  return probs;

}


// WINNER WINNER CHICKEN DINNER //
double exp_ggum_c(double theta,
                  double alpha,
                  double delta,
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

  double exp_1 = exp(alpha * (exp_mult * (theta - delta) - tau));

  return exp_1;

}


// [[Rcpp::export]]
NumericVector q_ggum(SEXP thetas,
                     SEXP params) {

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

  for(int p = 0; p < n_persons; p++){
    theta    = pthetas[p];
    exp_03   = exp_ggum_c(theta, alpha, delta, tau, 3);
    exp_12   = exp_ggum_c(theta, alpha, delta, tau, 1) +
               exp_ggum_c(theta, alpha, delta, tau, 2);
    probs[p] = (1 + exp_03) / (1 + exp_03 + exp_12);
  }

  UNPROTECT(2);

  return probs;

}

/*** R
library(microbenchmark)

thetas <- rnorm(100000)
params <- c(1, 0, -1)

microbenchmark(a <- q_ggum0(thetas, params),
               b <- q_ggum1(thetas, params),
               d <- q_ggum2(thetas, params),
               e <- q_ggum(thetas, params))

# q_ggum is faster
*/

////////////////
// Q GGUM ALL //
////////////////

// [[Rcpp::export]]
NumericMatrix q_ggum_all0(NumericMatrix thetas,
                          NumericMatrix params) {

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
    theta         = thetas(_, dim);
    param         = params(dim, _);
    probs(_, dim) = q_ggum(theta, param);
  }

  return probs;
}


double q_ggum_c(double theta,
                       double alpha,
                       double delta,
                       double tau) {

  // storage vectors
  double prob, exp_03, exp_12;

  // calculate parts of probability
  exp_03 = exp_ggum_c(theta, alpha, delta, tau, 3);
  exp_12 = exp_ggum_c(theta, alpha, delta, tau, 1) +
           exp_ggum_c(theta, alpha, delta, tau, 2);

  // put parts together to form probability
  prob   = (1 + exp_03) / (1 + exp_03 + exp_12);

  return prob;

}

// [[Rcpp::export]]
NumericMatrix q_ggum_all(NumericMatrix thetas,
                         NumericMatrix params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // dimensions of objects
  int n_persons = thetas.nrow();
  int n_items   = params.nrow();

  //declaring return
  NumericMatrix probs(Rf_allocMatrix(REALSXP, n_persons, n_items));

  for(int pers = 0; pers < n_persons; pers++){
    for(int item = 0; item < n_items; item++){
      probs[pers + item * n_persons] = q_ggum_c(thetas[pers + item * n_persons],
                                                params[item],
                                                params[item + 1 * n_items],
                                                params[item + 2 * n_items]);
    }
  }

  return probs;
}

/*** R
library(microbenchmark)

thetas <- cbind(rnorm(5000), rnorm(5000), rnorm(5000), rnorm(5000), rnorm(5000))
params <- rbind(c(1, 0, -1),
                c(.8, 2, 2),
                c(1.2, -1, 3),
                c(1, 1, 1),
                c(4, 3, 3))

microbenchmark(a <- q_ggum_all0(thetas, params),
               b <- q_ggum_all(thetas, params),
               d <- mupp:::q_ggum_all(thetas, params))

# q_ggum is faster
*/

///////////////
// QDER GGUM //
///////////////

// NumericVector pder1_theta_ggum(NumericVector thetas,
//                                NumericVector params) {
//
//   // Arguments:
//   //  - thetas: a vector of thetas across all people
//   //  - params: a vector of params for one item [alpha, delta, tau]
//   // Value:
//   //  - P'(Z_s = 1 | thetas)
//
//   // declaring parameters
//   double alpha = params[0];
//
//   // calculate exponent expressions
//   NumericVector exp_1 = exp_ggum(thetas, params, 1);
//   NumericVector exp_2 = exp_ggum(thetas, params, 2);
//   NumericVector exp_3 = exp_ggum(thetas, params, 3);
//
//   // calculate derivative of p1 (derivative of p0 is -dp1)
//   NumericVector dprobs = alpha * (exp_1 * (1 - 2 * exp_3) +
//                                   exp_2 * (2 - 1 * exp_3)) /
//                          pow((1 + exp_1 + exp_2 + exp_3), 2);
//
//   return dprobs;
// }
//
// // [[Rcpp::export]]
// NumericMatrix pder1_theta_ggum_all(NumericMatrix thetas,
//                                    NumericMatrix params) {
//
//   // Arguments:
//   //  - thetas: a matrix of thetas across all people/dims
//   //  - params: a matrix of params for one item [alpha, delta, tau]
//   // Value:
//   //  - P'(Z = 1 | thetas) across all dimensions
//
//   // declare number of persons/dimensions
//   int n_persons = thetas.nrow();
//   int n_dims    = params.nrow();
//
//   // indicating return vector
//   NumericMatrix dprobs(n_persons, n_dims);
//
//   // calculating probabilities across all dimensions
//   for(int dim = 0; dim < n_dims; dim++){
//     dprobs(_, dim) = pder1_theta_ggum(thetas(_, dim), params(dim, _));
//   }
//
//   return dprobs;
// }
//

/////////////////
// P MUPP PICK0 /
/////////////////

// [[Rcpp::export]]
NumericVector p_mupp_pick00(NumericMatrix Q,
                            int picked_dim = 1) {

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
NumericVector p_mupp_pick01(const NumericMatrix& Q,
                            int picked_dim = 1) {

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

# define get_dims(A) INTEGER(Rf_getAttrib(Q, R_DimSymbol));

// [[Rcpp::export]]
NumericVector p_mupp_pick0(const SEXP& Q,
                           int picked_dim = 1) {

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - picked_dim: the picked dimension
  // Value:
  //  - P(Z_s = 1 | thetas) individual element where s is the selected option

  // obtaining Q attributes
  int *dim_Q;

  dim_Q = get_dims(Q);

  // declare number of persons/dimensions
  int n_persons = 5000; //dim_Q[0];
  int n_dims    = 5; //dim_Q[1];

  // protect the objects
  PROTECT(Q);

  // pointers to objects
  double *pQ;
  pQ = REAL(Q);

  // vectors to store stuff
  NumericVector probs(n_persons, 1.0);
  double q;

  // calculating (intersection) probability across picked dimension
  for(int pers = 0; pers < n_persons; pers++){
    for(int dim = 0; dim < n_dims; dim++){
      q = pQ[pers + n_persons * dim];
      if(dim == picked_dim){
        probs[pers] *= (1 - q);
      } else{
        probs[pers] *= q;
      }
    }
  }

  UNPROTECT(1);

  return probs;
}

/*** R
library(microbenchmark)

thetas <- cbind(rnorm(5000), rnorm(5000), rnorm(5000), rnorm(5000), rnorm(5000))
params <- rbind(c(1, 0, -1),
                c(.8, 2, 2),
                c(1.2, -1, 3),
                c(1, 1, 1),
                c(4, 3, 3))
Q      <- q_ggum_all(thetas, params)

microbenchmark(a <- p_mupp_pick00(Q, 1),
               b <- p_mupp_pick01(Q, 1),
               d <- p_mupp_pick0(Q, 1))

# q_ggum is faster
*/

//
// NumericVector p_mupp_pick1(NumericMatrix Q,
//                            int picked_dim = 1) {
//
//   // Arguments:
//   //  - Q: a vector of probability of NOT selecting each option
//   //  - picked_dim: the picked dimension
//   // Value:
//   //  - P(Z_s = 1 | thetas) OVERALL where s is the selected option
//
//   // declare number of persons/dimensions
//   int n_persons = Q.nrow();
//   int n_dims    = Q.ncol();
//
//   // vectors to store stuff
//   NumericVector numer(n_persons);
//   NumericVector denom(n_persons);
//   NumericVector p(n_persons);
//
//   // calculating probability across picked dimension
//   for(int dim = 0; dim < n_dims; dim++){
//
//     // determine individual probability multiplication
//     p = p_mupp_pick0(Q, dim);
//
//     // add to numerator IF picked_dim, otherwise add to denominator
//     if(dim == picked_dim){
//       numer = numer + p;
//     } else{
//       denom = denom + p;
//     }
//   }
//
//   return numer / (numer + denom);
// }
//
// NumericVector p_mupp_rank1(NumericMatrix Q,
//                            IntegerVector order){
//
//   // Arguments:
//   //  - Q: a vector of probability of NOT selecting each option
//   //  - order: the order of the output
//   // Value:
//   //  - P(Z_s = 1 | thetas) where s is the selected option
//
//   // declare number of persons/dimensions
//   int n_persons = Q.nrow();
//   int n_dims    = order.size();
//
//   // matrices to store stuff
//   NumericMatrix Q_order = select_cols(Q, order);
//
//   // indicating return vector
//   NumericVector probs(n_persons, 1.0);
//
//   // calculating probability of being in particular order (1 by 1)
//   for(int dim = 0; dim < n_dims - 1; dim++){
//     probs = probs * p_mupp_pick1(Q_order(_, Range(dim, n_dims - 1)), 0);
//   }
//
//   return probs;
// }
//
// // [[Rcpp::export]]
// NumericMatrix p_mupp_rank_impl(NumericMatrix thetas,
//                                NumericMatrix params,
//                                IntegerVector dims            = NA_INTEGER,
//                                IntegerVector picked_order_id = NA_INTEGER) {
//
//   // Arguments:
//   //  - thetas: matrix of persons x dims (for all dims)
//   //  - params: matrix of dims x params (for single item dims)
//   //  - dims:   vector of dims of the items
//   //  - picked_order_id: index of picked order for each person
//   //                     (if NA, return ALL orders)
//   // Value:
//   //  - matrix of P(s > t > ...)(theta_s, theta_t, theta_ ...) for all s, t, ... in dims
//   //    (all of the possible permutations) OR
//   //  - matrix of P(s > t > ...)(theta_s, theta_t, theta_ ...) for the PICKED ORDER
//   //    (ONE permutation for each person)
//
//   // declare number of items/params (1) AND whether to return ALL combinations
//   int n_persons   = thetas.nrow();
//   int n_dims_all  = thetas.ncol();
//   int n_params    = params.nrow();
//   bool all_combs  = true;
//
//   // temporary vectors to store things
//   int n_orders;
//
//   // fixing dims if it is NA or missing
//   if(any(is_na(dims))){
//     dims = seq_len(n_params);
//   }
//
//   if(any(is_na(picked_order_id))){
//     picked_order_id = NA_INTEGER;
//   } else{
//     picked_order_id = rep_len(picked_order_id, n_persons);
//     all_combs       = false;
//   }
//
//   // R -> C Conversion (subtract 1)
//   IntegerVector dims_c            = dims - 1;
//   IntegerVector picked_order_id_c = picked_order_id - 1;
//
//   // declare number of items/params (2)
//   int n_dims_item = dims_c.size();
//
//   // Argument Checks (1):
//
//   // - dims must have the same number of elements as params
//   if(n_dims_item != n_params){
//     stop("dims must have the same number of elements as params");
//   }
//
//   // - dims must have MAX value less than thetas
//   if(max(dims_c) >= n_dims_all){
//     stop("dims must have max value less than the number of columns of theta");
//   }
//
//   // indicate temporary storage vectors
//   IntegerMatrix picked_orders = find_all_permutations(n_dims_item);
//   NumericMatrix Q             = q_ggum_all(select_cols(thetas, dims_c), params);
//
//   // Argument Checks (2):
//
//   // - picked_order_id must have MAX value at most factorial(n_dims) and MIN value at least 0
//   if(!all_combs){
//     if(max(picked_order_id_c) >= picked_orders.nrow()){
//       stop("picked_order_id must be at most the factorial(number of dimensions)");
//     } else if(min(picked_order_id_c) < 0){
//       stop("picked_order_id must be at least 1");
//     }
//   }
//
//   // determine the number of orders AND the Q matrix
//   //  - if we haven't specified selected order, we do this for EVERY ORDER
//   //  - otherwise, we do this JUST for the selected order
//   if(all_combs){
//     n_orders = picked_orders.nrow();
//   } else{
//     n_orders = 1;
//
//     // temporary vectors to store the Q_person and picked_order
//     NumericVector Q_person(n_params);
//     IntegerVector picked_order(n_params);
//
//     // for each person, rearrange Q so that [1, 2, 3, ...] is PICKED order ...
//     for(int person = 0; person < n_persons; person++){
//       picked_order = picked_orders(picked_order_id_c[person], _);
//       Q_person     = Q(person, _);
//       Q(person, _) = as<NumericVector>(Q_person[picked_order]);
//     }
//   }
//
//   // indicate return matrix
//   NumericMatrix probs(n_persons, n_orders);
//
//   // add all probabilities to return matrix
//   for(int perm = 0; perm < n_orders; perm++){
//     probs(_, perm) = p_mupp_rank1(Q, picked_orders(perm, _));
//   }
//
//   return probs;
// }
//
// //[[Rcpp::export]]
// NumericMatrix loglik_mupp_rank_impl(NumericMatrix thetas,
//                                     NumericMatrix params,
//                                     IntegerMatrix items,
//                                     IntegerMatrix picked_orders){
//
//   // Arguments:
//   //  - thetas: matrix of persons x dims (for all dims)
//   //  - params: matrix of dims x params (for all params, in order)
//   //  - items: matrix of [item, statement, dim] for all items, where
//   //           statement is integer indicating the statement number and aligns
//   //           with params
//   //  - picked_order_id: matrix if picked_order_id for [people x items]
//
//   // Value:
//   //  - matrix of loglikelihoods for mupp rank stuff ... (potentially aggregated)
//
//   // declare number of persons and dimensions of everything else
//   int n_persons    = thetas.nrow();
//   int n_dims_all   = thetas.ncol();
//   int n_statements = params.nrow();
//   int n_params     = params.ncol();
//   int n_rows_items = items.nrow();
//   int n_cols_items = items.ncol();
//
//   int n_rows_resp  = picked_orders.nrow();
//   int n_cols_resp  = picked_orders.ncol();
//
//   // Argument Checks (1)
//
//   // - items must have three columns
//   if(n_cols_items != 3){
//     stop("items must have three columns [item, statement, dimension]");
//   }
//
//   // - params must have three columns
//   if(n_params != 3){
//     stop("params must have three columns [alpha, delta, tau]");
//   }
//
//   // - number of response rows must match number of persons
//   if(n_persons != n_rows_resp){
//     stop("number of persons must match number of rows of the resp matrix");
//   }
//
//   // declare items
//   IntegerVector item_ids      = items(_, 0);
//   IntegerVector statement_ids = items(_, 1) - 1;
//   IntegerVector dim_ids       = items(_, 2);
//   IntegerVector unique_items  = sort_unique(item_ids);
//   int n_items = unique_items.size();
//
//   // Argument Checks (2)
//
//   // - number of response cols must match number of items
//   if(n_items != n_cols_resp){
//     stop("number of unique item ids must match number of cols of the resp matrix");
//   }
//
//   // - dims must be between 1 and number of implied dims
//   if((min(dim_ids) <= 0) | (max(dim_ids) > n_dims_all)){
//     stop("number of dimensions in item matrix must be between 1 and the number of columns of theta");
//   }
//
//   // - statements must be between 1 and number of implied statements
//   if((min(statement_ids) < 0) | (max(statement_ids) >= n_statements)){
//     stop("number of statements in item matrix must be between 1 and the number or rows of params");
//   }
//
//   // indicate return matrix and temporary vectors
//   NumericMatrix loglik(n_persons, n_items);
//   NumericVector p(n_persons);
//   LogicalVector item_flag(n_rows_items);
//
//   for(int item = 0; item < n_items; item++){
//
//     // pulling out statements and dimensions for this particular item
//     item_flag = (item_ids == unique_items[item]);
//
//     // calculating probability (all thetas, relevant params/dims, item orders)
//     p         = p_mupp_rank_impl(thetas,
//                                  select_rows(params, statement_ids[item_flag]),
//                                  dim_ids[item_flag],
//                                  picked_orders(_, item))(_, 0);
//
//   // calculating likelihood and putting it in matrix
//     loglik(_, item) = log(p);
//   }
//
//   return loglik;
// }
