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
#include "ggum_probs.h"
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


// arrange matrix so that every row is sorted by picked order
void arrange_by_picked(NumericMatrix M,
                       const IntegerVector & indices,
                       const IntegerMatrix & orders,
                       bool reverse = false) {

  // declare number of rows/columns
  int n_cols = M.ncol();
  int n_rows = M.nrow();

  // temporary vectors to store the matrices and picked order
  NumericVector M_person     = no_init(n_cols);
  IntegerVector picked_order = no_init(n_cols);
  IntegerMatrix all_orders   = clone(orders);

  // indices to reverse the ordering (by row) if required)
  if(reverse){
    for(int row = 0; row < all_orders.nrow(); row++){
      IntegerVector row_order   = all_orders.row(row);
      IntegerVector match_order = clone(row_order).sort();
      all_orders.row(row)       = match(match_order, row_order) - 1;
    }
  }

  // for each person, rearrange M so that [1, 2, 3] is the PICKED order ...
  for(int person = 0; person < n_rows; person++){

    // pull out appropriate order
    picked_order  = all_orders.row(indices[person]);

    // rearrange row and assign to temporary vector
    M_person      = M.row(person);

    // put back in place of matrix
    M.row(person) = as<NumericVector>(M_person[picked_order]);
  }

}

// saved permutations (so we don't have to recreate this all of the time)
List saved_permutations = List::create(find_all_permutations(1),
                                       find_all_permutations(2),
                                       find_all_permutations(3),
                                       find_all_permutations(4),
                                       find_all_permutations(5),
                                       find_all_permutations(6));

// extract save permutation (if we have it stored)
// [[Rcpp::export]]
IntegerMatrix extract_permutations(int n,
                                   int init = 0){
  if(n < saved_permutations.size() & init == 0){
    return saved_permutations[n - 1];
  } else{
    return find_all_permutations(n, init);
  }
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
      probs = probs * (1 - Q.column(dim));
    } else{
      probs = probs * Q.column(dim);
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
NumericMatrix pder1_mupp_rank1(const NumericMatrix & P,
                               const NumericMatrix & dP,
                               const IntegerVector order){

  // Arguments:
  //  - P:  a matrix of probability of selecting each option
  //  - dP: derivative of the probability matrix (same dimensions as P)
  //  - order: the order of the output
  // Value:
  //  - P'(Z_s = 1 | thetas) for s and t

  // At some point - create pder1_mupp_pick0, pder1_mupp_pick1, ...
  //               - should be doable with 3+ dimensions, as pick is relatively
  //                 simple, so the derivative should be simple, and rank is just
  //                 pick multiplied: A, B, C --> dP/dt1 = dA * BC (t1 isn't in B/C)
  //                                          --> dP/dt2 = d(AB) * C, ...

  // declare number of persons/dimensions
  int n_persons  = P.nrow();
  int n_dims_all = P.ncol();
  int n_dims     = order.size();

  // argument checks
  if(n_persons != dP.nrow()){
    stop("P and dP do not have the same number of rows");
  }

  if(n_dims_all != dP.ncol()){
    stop("P and dP do not have the same number of columns");
  }

  if(n_dims != 2){
    stop("pder1_mupp_rank1 only provides derivatives for 2 dimensions, so far");
  }

  // indicating matrices to store things
  NumericMatrix P_o  = select_cols(P, order);
  NumericMatrix dP_o = select_cols(dP, order);
  int dim, alt, out;

  // indicating return matrix
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // fixing probabilities/derivatives for NOT chosen categories
  P_o.column(1)  = 1 - P_o.column(1);
  dP_o.column(1) = -dP_o.column(1);

  // denominator of the model
  NumericVector denom = P_o.column(0) * P_o.column(1) + (1 - P_o.column(0)) * (1 - P_o.column(1));

  // set dprobs to appropriate column of derivative matrix
  for(dim = 0; dim < n_dims; dim++){
    alt                = n_dims - dim - 1;
    out                = order[dim];
    dprobs.column(out) = P_o.column(alt) * (1 - P_o.column(alt)) * dP_o.column(dim) / pow(denom, 2);
  }

  return dprobs;

}

List initialize_mupp(const NumericMatrix & thetas,
                     const NumericMatrix & params,
                     IntegerVector dims            = NA_INTEGER,
                     IntegerVector picked_order_id = NA_INTEGER) {

  // declare number of items/params (1) AND whether to return ALL combinations
  int n_persons   = thetas.nrow();
  int n_dims_all  = thetas.ncol();
  int n_params    = params.nrow();
  bool all_combs  = true;

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

  // Argument Checks (2):

  // - picked_order_id must have MAX value at most factorial(n_dims) and MIN value at least 0
  if(!all_combs){
    if(max(picked_order_id_c) >= picked_orders.nrow()){
      stop("picked_order_id must be at most the factorial(number of dimensions)");
    } else if(min(picked_order_id_c) < 0){
      stop("picked_order_id must be at least 1");
    }
  }

  return List::create(
    _["dims"]            = dims_c,
    _["picked_order_id"] = picked_order_id_c,
    _["picked_orders"]   = picked_orders,
    _["n_persons"]       = n_persons,
    _["n_params"]        = n_params,
    _["n_dims_all"]      = n_dims_all,
    _["n_dims_item"]     = n_dims_item,
    _["all_combs"]       = all_combs
  );
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

  // initialize parameters and return stuff
  List all_params  = initialize_mupp(thetas, params, dims, picked_order_id);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"];
  bool all_combs                  = all_params["all_combs"];
  IntegerVector dims_c            = all_params["dims"];
  IntegerVector picked_order_id_c = all_params["picked_order_id"];
  IntegerMatrix picked_orders     = all_params["picked_orders"];

  // indicate temporary storage vectors
  int n_orders;
  NumericMatrix Q  = q_ggum_all(select_cols(thetas, dims_c), params);

  // determine the number of orders AND the Q matrix
  //  - if we haven't specified selected order, we do this for EVERY ORDER
  //  - otherwise, we do this JUST for the selected order
  if(all_combs){
    n_orders = picked_orders.nrow();
  } else{
    n_orders = 1;
    arrange_by_picked(Q,
                      picked_order_id_c,
                      picked_orders);
  }

  // indicate return matrix
  NumericMatrix probs(n_persons, n_orders);

  // add all probabilities to return matrix
  for(int perm = 0; perm < n_orders; perm++){
    probs.column(perm) = p_mupp_rank1(Q, picked_orders.row(perm));
  }

  return probs;
}


// [[Rcpp::export]]
List pder1_mupp_rank_impl(NumericMatrix & thetas,
                          const NumericMatrix & params,
                          IntegerVector dims            = NA_INTEGER,
                          IntegerVector picked_order_id = NA_INTEGER) {

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for single item dims)
  //  - dims:   vector of dims of the items
  //  - picked_order_id: index of picked order for each person
  //                     (if NA, return ALL orders)
  //
  // Value:
  //  - list of P'(s > t > ...)(theta_s, theta_t, theta_ ...) ds, where the list
  //    elements correspond to the choice and the column corresponds to the
  //    dimension

  // initialize parameters and return stuff
  List all_params  = initialize_mupp(thetas, params, dims, picked_order_id);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"];
  int n_dims_all                  = all_params["n_dims_all"];
  int n_dims_item                 = all_params["n_dims_item"];
  bool all_combs                  = all_params["all_combs"];
  IntegerVector dims_c            = all_params["dims"];
  IntegerVector picked_order_id_c = all_params["picked_order_id"];
  IntegerMatrix picked_orders     = all_params["picked_orders"];

  // update thetas to be in the correct order and just those selected
  thetas           = select_cols(thetas, dims_c);

  // indicate temporary storage vectors
  int n_orders;
  NumericMatrix P  = p_ggum_all(thetas, params);
  NumericMatrix dP = pder1_ggum_all(thetas, params);

  // determine the number of orders AND the Q matrix
  //  - if we haven't specified selected order, we do this for EVERY ORDER
  //  - otherwise, we do this JUST for the selected order
  if(all_combs){
    n_orders = picked_orders.nrow();
  } else{
    n_orders = 1;
    arrange_by_picked(P,
                      picked_order_id_c,
                      picked_orders);
    arrange_by_picked(dP,
                      picked_order_id_c,
                      picked_orders);
  }

  // indicate return list
  List dprobs_all;

  // add all probabilities to return matrix
  for(int perm = 0; perm < n_orders; perm++){

    // initialize probability matrix to 0 and calculate probabilities
    NumericMatrix dprobs(n_persons, n_dims_all);
    NumericMatrix dprobs_ = pder1_mupp_rank1(P, dP, picked_orders.row(perm));

    // reverse ordering to get back to appropriate order (given dims)
    if(!all_combs){
      arrange_by_picked(dprobs_,
                        picked_order_id_c,
                        picked_orders,
                        true);
    }

    // add probabilities to appropriate column of matrix (combining same dims)
    for(int dim = 0; dim < n_dims_item; dim++){
      dprobs.column(dims_c[dim]) = dprobs.column(dims_c[dim]) + dprobs_.column(dim);
    }

    // assign back to list
    dprobs_all.push_back(dprobs);

  }

  return dprobs_all;
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
  IntegerVector item_ids      = items.column(0);
  IntegerVector statement_ids = items.column(1) - 1;
  IntegerVector dim_ids       = items.column(2);
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
                                 picked_orders.column(item)).column(0);

  // calculating likelihood and putting it in matrix
    loglik.column(item) = log(p);
  }

  return loglik;
}
