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
#include "mupp_utilities.h"
using namespace Rcpp;

/* HELPER FUNCTIONS
 *   - select_cols (select non-proximate columns)
 *   - select_rows (select non-proximate rows)
 *   - extract_permutations
 */

NumericMatrix select_cols(const NumericMatrix & X,
                          const IntegerVector & ind) {

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

// saved permutations (so we don't have to recreate this all of the time)
List saved_permutations = List::create(find_all_permutations(1),
                                       find_all_permutations(2),
                                       find_all_permutations(3),
                                       find_all_permutations(4),
                                       find_all_permutations(5),
                                       find_all_permutations(6));

IntegerMatrix extract_permutations(int n,
                                   int init = 0){

  // pull out saved_permutation OR calculate (if not existing)
  if(n < saved_permutations.size() & init == 0){
    return saved_permutations[n - 1];
  } else{
    return find_all_permutations(n, init);
  }

}

NumericMatrix select_rows(const NumericMatrix & X,
                          const IntegerVector & ind) {

  // Arguments:
  //  - X:   a matrix to select rows
  //  - ind: an integer vector of rows to select

  // capturing constants
  int out_rows = ind.size();
  int out_cols = X.ncol();
  int inp_rows = X.nrow();

  NumericMatrix Y(Rf_allocMatrix(REALSXP, out_rows, out_cols));

  for(int j = 0; j < out_cols; j++){
    for(int i = 0; i < out_rows; i++){
      Y[i + j * out_rows] = X[ind[i] + j * inp_rows];
    }
  }

  return Y;
}

/* GGUM STUFF
 *  - exp_ggum/exp_ggum_c
 *  - p_ggum/q_ggum_all
 *  - p_der1_theta_ggum
 */

NumericVector exp_ggum(const SEXP & thetas,
                       const SEXP & params,
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

double exp_ggum_c(const double theta,
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

  double exp_1 = exp(alpha * (exp_mult * (theta - delta) - tau));

  return exp_1;

}

NumericVector q_ggum(const SEXP & thetas,
                     const SEXP & params) {

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
    exp_03   = exp_ggum_c(theta, alpha, delta, tau, 3);
    exp_12   = exp_ggum_c(theta, alpha, delta, tau, 1) +
               exp_ggum_c(theta, alpha, delta, tau, 2);
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
    theta         = thetas(_, dim);
    param         = params(dim, _);
    probs(_, dim) = q_ggum(theta, param);
  }

  return probs;
}

NumericVector pder1_theta_ggum(const NumericVector & thetas,
                               const NumericVector & params) {

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

// [[Rcpp::export]]
NumericMatrix pder1_theta_ggum_all(const NumericMatrix & thetas,
                                   const NumericMatrix & params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P'(Z = 1 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow();
  int n_dims    = params.nrow();

  // indicating return vector
  NumericMatrix dprobs(n_persons, n_dims);

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    dprobs(_, dim) = pder1_theta_ggum(thetas(_, dim), params(dim, _));
  }

  return dprobs;
}

/* MUPP STUFF
 *   - p_mupp_pick0     - numerator (individual component) of PICK probability
 *   - p_mupp_pick1     - OVERALL PICK probability
 *   - p_mupp_rank1     - OVERALL RANK probability
 *   - p_mupp_rank_impl - OVERALL RANK probability for an individual order OR
 *                        all possible order combinations
 */

NumericVector p_mupp_pick0(const NumericMatrix & Q,
                           const int picked_dim = 1) {

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

NumericVector p_mupp_pick1(const NumericMatrix & Q,
                           const int picked_dim = 1) {

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - picked_dim: the picked dimension
  // Value:
  //  - P(Z_s = 1 | thetas) OVERALL where s is the selected option

  // declare number of persons/dimensions
  int n_persons = Q.nrow();
  int n_dims    = Q.ncol();

  // vectors to store stuff
  NumericVector numer(n_persons, 0.0);
  NumericVector denom(n_persons, 0.0);
  NumericVector p = no_init(n_persons);

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

NumericVector p_mupp_rank1(const NumericMatrix & Q,
                           const IntegerVector order){

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

// [[Rcpp::export]]
NumericMatrix p_mupp_rank_impl(const NumericMatrix & thetas,
                               const NumericMatrix & params,
                               IntegerVector dims            = NA_INTEGER,
                               IntegerVector picked_order_id = NA_INTEGER) {

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for single item dims)
  //  - dims:   vector of dims of the items
  //  - picked_order_id: index of picked order for each person
  //                     (if NA, return ALL orders)
  // Value:
  //  - matrix of P(s > t > ...)(theta_s, theta_t, theta_ ...) for all s, t, ... in dims
  //    (all of the possible permutations) OR
  //  - matrix of P(s > t > ...)(theta_s, theta_t, theta_ ...) for the PICKED ORDER
  //    (ONE permutation for each person)

  // declare number of items/params (1) AND whether to return ALL combinations
  int n_persons   = thetas.nrow();
  int n_dims_all  = thetas.ncol();
  int n_params    = params.nrow();
  bool all_combs  = true;

  // temporary vectors to store things
  int n_orders;

  // forcing dims and picked order ID to have positive length
  if(dims.size() == 0){
    dims = IntegerVector::create(NA_INTEGER);
  }

  if(picked_order_id.size() == 0){
    picked_order_id = IntegerVector::create(NA_INTEGER);
  }

  // fixing dims if it is NA or missing
  if(dims[0] == NA_INTEGER){
    dims = seq_len(n_params);
  }

  if(picked_order_id[0] == NA_INTEGER){
    picked_order_id = NA_INTEGER;
  } else{
    picked_order_id = rep_len(picked_order_id, n_persons);
    all_combs       = false;
  }

  // R -> C Conversion (subtract 1)
  IntegerVector dims_c            = dims - 1;
  IntegerVector picked_order_id_c = picked_order_id - 1;

  // declare number of items/params (2)
  int n_dims_item = dims_c.size();

  // Argument Checks (1):

  // - dims must have the same number of elements as params
  if(n_dims_item != n_params){
    stop("dims must have the same number of elements as params");
  }

  // - dims must have MAX value less than thetas
  if(max(dims_c) >= n_dims_all){
    stop("dims must have max value less than the number of columns of theta");
  }

  // indicate temporary storage vectors
  IntegerMatrix picked_orders = extract_permutations(n_dims_item);
  NumericMatrix Q             = q_ggum_all(select_cols(thetas, dims_c), params);

  // Argument Checks (2):

  // - picked_order_id must have MAX value at most factorial(n_dims) and MIN value at least 0
  if(!all_combs){
    if(max(picked_order_id_c) >= picked_orders.nrow()){
      stop("picked_order_id must be at most the factorial(number of dimensions)");
    } else if(min(picked_order_id_c) < 0){
      stop("picked_order_id must be at least 1");
    }
  }

  // determine the number of orders AND the Q matrix
  //  - if we haven't specified selected order, we do this for EVERY ORDER
  //  - otherwise, we do this JUST for the selected order
  if(all_combs){
    n_orders = picked_orders.nrow();
  } else{
    n_orders = 1;

    // temporary vectors to store the Q_person and picked_order
    NumericVector Q_person     = no_init(n_params);
    IntegerVector picked_order = no_init(n_params);

    // for each person, rearrange Q so that [1, 2, 3, ...] is PICKED order ...
    for(int person = 0; person < n_persons; person++){
      picked_order = picked_orders(picked_order_id_c[person], _);
      Q_person     = Q(person, _);
      Q(person, _) = as<NumericVector>(Q_person[picked_order]);
    }
  }

  // indicate return matrix
  NumericMatrix probs(n_persons, n_orders);

  // add all probabilities to return matrix
  for(int perm = 0; perm < n_orders; perm++){
    probs(_, perm) = p_mupp_rank1(Q, picked_orders(perm, _));
  }

  return probs;
}

//[[Rcpp::export]]
NumericMatrix loglik_mupp_rank_impl(const NumericMatrix & thetas,
                                    const NumericMatrix & params,
                                    const IntegerMatrix & items,
                                    const IntegerMatrix & picked_orders){

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for all params, in order)
  //  - items: matrix of [item, statement, dim] for all items, where
  //           statement is integer indicating the statement number and aligns
  //           with params
  //  - picked_order_id: matrix if picked_order_id for [people x items]

  // Value:
  //  - matrix of loglikelihoods for mupp rank stuff ... (potentially aggregated)

  // declare number of persons and dimensions of everything else
  int n_persons    = thetas.nrow();
  int n_dims_all   = thetas.ncol();
  int n_statements = params.nrow();
  int n_params     = params.ncol();
  int n_rows_items = items.nrow();
  int n_cols_items = items.ncol();

  int n_rows_resp  = picked_orders.nrow();
  int n_cols_resp  = picked_orders.ncol();

  // Argument Checks (1)

  // - items must have three columns
  if(n_cols_items != 3){
    stop("items must have three columns [item, statement, dimension]");
  }

  // - params must have three columns
  if(n_params != 3){
    stop("params must have three columns [alpha, delta, tau]");
  }

  // - number of response rows must match number of persons
  if(n_persons != n_rows_resp){
    stop("number of persons must match number of rows of the resp matrix");
  }

  // declare items
  IntegerVector item_ids      = items(_, 0);
  IntegerVector statement_ids = items(_, 1) - 1;
  IntegerVector dim_ids       = items(_, 2);
  IntegerVector unique_items  = sort_unique(item_ids);
  int n_items = unique_items.size();

  // Argument Checks (2)

  // - number of response cols must match number of items
  if(n_items != n_cols_resp){
    stop("number of unique item ids must match number of cols of the resp matrix");
  }

  // - dims must be between 1 and number of implied dims
  if((min(dim_ids) <= 0) | (max(dim_ids) > n_dims_all)){
    stop("number of dimensions in item matrix must be between 1 and the number of columns of theta");
  }

  // - statements must be between 1 and number of implied statements
  if((min(statement_ids) < 0) | (max(statement_ids) >= n_statements)){
    stop("number of statements in item matrix must be between 1 and the number or rows of params");
  }

  // indicate return matrix and temporary vectors
  NumericMatrix loglik(Rf_allocMatrix(REALSXP, n_persons, n_items));
  NumericVector p         = no_init(n_persons);
  LogicalVector item_flag = no_init(n_rows_items);

  for(int item = 0; item < n_items; item++){

    // pulling out statements and dimensions for this particular item
    item_flag = (item_ids == unique_items[item]);

    // calculating probability (all thetas, relevant params/dims, item orders)
    p         = p_mupp_rank_impl(thetas,
                                 select_rows(params, statement_ids[item_flag]),
                                 dim_ids[item_flag],
                                 picked_orders(_, item))(_, 0);

  // calculating likelihood and putting it in matrix
    loglik(_, item) = log(p);
  }

  return loglik;
}
