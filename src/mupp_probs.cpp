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

/* MUPP PROBABILITY STUFF
 *   - p_mupp_pick0     - numerator (individual component) of PICK probability
 *   - p_mupp_pick1     - OVERALL PICK probability
 *   - p_mupp_rank1     - OVERALL RANK probability
 *   - p_mupp_rank_impl - OVERALL RANK probability for an individual order OR
 *                        all possible order combinations
 *   - pder1_mupp_rank_impl - DERIVATIVES for an individual order OR all possible
 *                            combinations
 */
NumericVector p_mupp_pick0(const NumericMatrix & Q,
                           const int picked_dim = 1) {

  // Arguments:
  //  - Q: a vector of probability of NOT selecting each option
  //  - picked_dim: the picked dimension
  // Value:
  //  - P(Z_s = 1 | thetas) individual element where s is the selected option

  // declare number of persons/dimensions
  int n_persons = Q.nrow(),
      n_dims    = Q.ncol();

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

  // protection for really rare events
  LogicalVector flags = (probs == 0) | (probs == 1);

  if(is_true(any(flags))){
    probs[probs == 0] = 0 + 1e-09;
    probs[probs == 1] = 1 - 1e-09;
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
  int n_persons = Q.nrow(),
      n_dims    = Q.ncol();

  // vectors to store stuff
  NumericVector numer(n_persons, 0.0),
                denom(n_persons, 0.0),
                p = no_init(n_persons);

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
  int n_persons = Q.nrow(),
      n_dims    = order.size();

  // matrices to store stuff
  NumericMatrix Q_order = select_cols(Q, order);

  // indicating return vector
  NumericVector probs(n_persons, 1.0);

  // calculating probability of being in particular order (1 by 1)
  for(int dim = 0; dim < n_dims - 1; dim++){
    probs = probs * p_mupp_pick1(Q_order(_, Range(dim, n_dims - 1)), 0);
  }

  // protection for really rare events
  LogicalVector flags = (probs == 0) | (probs == 1);

  if(is_true(any(flags))){
    probs[probs == 0] = 0 + 1e-09;
    probs[probs == 1] = 1 - 1e-09;
  }

  return probs;
}

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
  int n_persons  = P.nrow(),
      n_dims_all = P.ncol(),
      n_dims     = order.size();

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
  NumericMatrix P_o  = select_cols(P, order),
                dP_o = select_cols(dP, order);
  int dim, alt, out;

  // indicating return matrix
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims));

  // fixing probabilities/derivatives for NOT chosen categories
  P_o.column(1)  = -P_o.column(1) + 1;
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

NumericMatrix pder2_mupp_rank1(const NumericMatrix & P,
                               const NumericMatrix & dP,
                               const NumericMatrix & d2P,
                               const IntegerVector order){

  // Arguments:
  //  - P:     a matrix of probability of selecting each option
  //  - dP:    derivative of the probability matrix (same dimensions as P)
  //  - d2P:   second derivative of the probability matrix (same dimensions as P)
  //  - order: the order of the output
  // Value:
  //  - P''(Z_s = 1 | thetas) for s and t and st
  //    (generalize to more than two dimensions?)

  // declare number of persons/dimensions
  int n_persons  = P.nrow(),
      n_dims_all = P.ncol(),
      n_dims     = order.size(),
      n_pairs    = (n_dims * (n_dims - 1)) / 2;

  // argument checks
  if(n_persons != dP.nrow()){
    stop("P and dP do not have the same number of rows");
  }

  if(n_dims_all != dP.ncol()){
    stop("P and dP do not have the same number of columns");
  }

  if(n_dims != 2){
    stop("pder2_mupp_rank1 only provides derivatives for 2 dimensions, so far");
  }

  // indicating matrices to store things
  NumericMatrix P_o   = select_cols(P,   order),
                dP_o  = select_cols(dP,  order),
                d2P_o = select_cols(d2P, order);
  NumericVector p1   = no_init(n_persons),
                p2   = no_init(n_persons),
                dp1  = no_init(n_persons),
                dp2  = no_init(n_persons),
                d2p1 = no_init(n_persons);
  int dim, alt, out;

  // indicating return matrix
  NumericMatrix dprobs(Rf_allocMatrix(REALSXP, n_persons, n_dims + n_pairs));

  // fixing probabilities/derivatives for NOT chosen categories
  P_o.column(1)   = -P_o.column(1) + 1;
  dP_o.column(1)  = -dP_o.column(1);
  d2P_o.column(1) = -d2P_o.column(1);

  // denominator of the model
  NumericVector denom  = P_o.column(0) * P_o.column(1) + (1 - P_o.column(0)) * (1 - P_o.column(1));
  NumericVector denom3 = pow(denom, 3);

  // set dprobs to appropriate column of derivative matrix
  for(dim = 0; dim < n_dims; dim++){

    // determining dimensions
    alt   = n_dims - dim - 1;
    out   = order[dim];

    // assigning parts
    p1   = P_o.column(dim);  p2  = P_o.column(alt);
    dp1  = dP_o.column(dim); dp2 = dP_o.column(alt);
    d2p1 = d2P_o.column(dim);

    // calculating (second derivative)
    dprobs.column(out) = (p2 * (1 - p2)) *
                         (d2p1 * denom - 2 * pow(dp1, 2) * (2 * p2 - 1)) /
                         denom3;

    // calculating (cross derivative in third column)
    if(dim == 1){
      dprobs.column(dim + 1) = dp1 * dp2 * (1 - p1 - p2) / denom3;
    }
  }

  return dprobs;
}

List initialize_mupp_p(const NumericMatrix & thetas,
                       const NumericMatrix & params,
                       IntegerVector dims            = 0,
                       IntegerVector picked_order_id = 0) {

  // declare number of items/params (1) AND whether to return ALL combinations
  int n_persons  = thetas.nrow(),
      n_dims_all = thetas.ncol(),
      n_params   = params.nrow();
  bool all_combs = true;

  // forcing dims and picked order ID to have positive length
  if(dims.size() == 0){
    dims = IntegerVector::create(0);
  }

  if(picked_order_id.size() == 0){
    picked_order_id = IntegerVector::create(0);
  }

  // fixing dims if it is NA or missing
  if(dims[0] == 0){
    dims = seq_len(n_params);
  }

  if((picked_order_id[0] == 0) | is_true(all(is_na(picked_order_id)))){
    picked_order_id = NA_INTEGER;
  } else{
    picked_order_id = rep_len(picked_order_id, n_persons);
    all_combs       = false;
  }

  // R -> C Conversion (subtract 1)
  IntegerVector dims_c            = dims - 1,
                picked_order_id_c = picked_order_id - 1;

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
    } else if(min(na_omit(picked_order_id_c)) < 0){
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
                               IntegerVector dims            = 0,
                               IntegerVector picked_order_id = 0) {

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
  List all_params  = initialize_mupp_p(thetas, params, dims, picked_order_id);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"];
  bool all_combs                  = all_params["all_combs"];
  IntegerVector dims_c            = all_params["dims"],
                picked_order_id_c = all_params["picked_order_id"];
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
ListOf<NumericMatrix> pder1_mupp_rank_impl(const NumericMatrix & thetas,
                                           const NumericMatrix & params,
                                           IntegerVector dims            = 0,
                                           IntegerVector picked_order_id = 0) {

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for single item dims)
  //  - dims:   vector of dims of the items
  //  - picked_order_id: index of picked order for each person
  //                     (if 0, return ALL orders)
  //
  // Value:
  //  - list of P'(s > t > ...)(theta_s, theta_t, theta_ ...) ds, where the list
  //    elements correspond to the choice and the column corresponds to the
  //    dimension

  // initialize parameters and return stuff
  List all_params  = initialize_mupp_p(thetas, params, dims, picked_order_id);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"],
      n_dims_all                  = all_params["n_dims_all"],
      n_dims_item                 = all_params["n_dims_item"];
  bool all_combs                  = all_params["all_combs"];
  IntegerVector dims_c            = all_params["dims"],
                picked_order_id_c = all_params["picked_order_id"];
  IntegerMatrix picked_orders     = all_params["picked_orders"];

  // update thetas to be in the correct order and just those selected
  NumericMatrix thetas_ordered    = select_cols(thetas, dims_c);

  // indicate temporary storage vectors
  int n_orders;
  NumericMatrix P  = p_ggum_all(thetas_ordered, params),
                dP = pder1_ggum_all(thetas_ordered, params);

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
    NumericMatrix dprobs(n_persons, n_dims_all),
                  dprobs_ = pder1_mupp_rank1(P, dP, picked_orders.row(perm));

    // reverse ordering to get back to appropriate order (given dims)
    if(!all_combs){
      arrange_by_picked(dprobs_,
                        picked_order_id_c,
                        picked_orders,
                        true);
    }

    // add deriv to appropriate column of matrix (combining same dims)
    for(int dim = 0; dim < n_dims_item; dim++){
      dprobs.column(dims_c[dim]) = dprobs.column(dims_c[dim]) + dprobs_.column(dim);
    }

    // assign back to list
    dprobs_all.push_back(dprobs);

  }

  return dprobs_all;
}

// [[Rcpp::export]]
ListOf<NumericMatrix> pder2_mupp_rank_impl(const NumericMatrix & thetas,
                                           const NumericMatrix & params,
                                           IntegerVector dims            = 0,
                                           IntegerVector picked_order_id = 0) {

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for single item dims)
  //  - dims:   vector of dims of the items
  //  - picked_order_id: index of picked order for each person
  //                     (if NA, return ALL orders)
  //
  // Value:
  //  - list of P''(s > t > ...)(theta_s, theta_t, theta_ ...) ds, where the list
  //    elements correspond to the choice, the column corresponds to the
  //    dimension, and the remaining columns corresponds to the cross-dimension

  // initialize parameters and return stuff
  List all_params  = initialize_mupp_p(thetas, params, dims, picked_order_id);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"],
      n_dims_all                  = all_params["n_dims_all"],
      n_dims_item                 = all_params["n_dims_item"];
  bool all_combs                  = all_params["all_combs"];
  IntegerVector dims_c            = all_params["dims"],
                picked_order_id_c = all_params["picked_order_id"];
  IntegerMatrix picked_orders     = all_params["picked_orders"];

  // update thetas to be in the correct order and just those selected
  NumericMatrix thetas_ordered    = select_cols(thetas, dims_c);

  // number of pairs (off-diagonal)
  int n_pairs = n_dims_all * (n_dims_all - 1) / 2;
  int dim, dim1, dim2, dim12, cross_mult;

  // NOTE: // * // indicates things that need to be updated if ranking more than 2

  // indicate temporary storage vectors
  int n_orders;
  NumericMatrix P   = p_ggum_all(thetas_ordered, params),
                dP  = pder1_ggum_all(thetas_ordered, params),
                d2P = pder2_ggum_all(thetas_ordered, params);

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
    arrange_by_picked(d2P,
                      picked_order_id_c,
                      picked_orders);
  }

  // indicate return list
  List dprobs_all;

  // add all probabilities to return matrix
  for(int perm = 0; perm < n_orders; perm++){

    // initialize probability matrix to 0 and calculate probabilities
    NumericMatrix dprobs(n_persons, n_dims_all + n_pairs),
                  dprobs_ = pder2_mupp_rank1(P, dP, d2P, picked_orders.row(perm));

    // * // reverse ordering to get back to appropriate order (given dims)
    //      (need to reverse-order cross-products)
    if(!all_combs){
      arrange_by_picked(dprobs_,
                        picked_order_id_c,
                        picked_orders,
                        true);
    }

    // add second deriv to appropriate column of matrix (combining same dims)
    for(dim = 0; dim < n_dims_item; dim++){
      dim1 = dims_c[dim];
      dprobs.column(dim1) = dprobs.column(dim1) + dprobs_.column(dim);
    }

    // add cross-derivative of probability
    for(int d1 = 0; d1 < n_dims_item - 1; d1++){
      for(int d2 = 1; d2 < n_dims_item; d2++){

        // add 1 to dim (start at cross-prod elements)
        dim += 1;

        // add cross-derivative of probability
        dim1 = dims_c[d1];
        dim2 = dims_c[d2];

        // dim1 =  dim2, everything (both off diags) goes into elt
        // dim1 != dim2, include one of the off-diags
        if(dim1 == dim2){
          cross_mult = 2;
        } else{
          cross_mult = 1;
        }

        // column to put the cross product in
        dim12 = find_crossprod_column(dim1, dim2, n_dims_all);

        // add cross-deriv in appropriate column of matrix
        dprobs.column(dim12) = dprobs.column(dim12) + cross_mult * dprobs_.column(2);

        // assign back to list
        dprobs_all.push_back(dprobs);
      }
    }
  }

  return dprobs_all;
}


/* MUPP LIKELIHOOD STUFF
 */

List initialize_mupp_l(const NumericMatrix & thetas,
                       const NumericMatrix & params,
                       const IntegerMatrix & items,
                       const IntegerMatrix & picked_orders) {

  // declare number of persons and dimensions of everything else
  int n_persons    = thetas.nrow(),
      n_dims_all   = thetas.ncol(),
      n_statements = params.nrow(),
      n_params     = params.ncol(),
      n_rows_items = items.nrow(),
      n_cols_items = items.ncol();

  int n_rows_resp  = picked_orders.nrow(),
      n_cols_resp  = picked_orders.ncol();

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
  IntegerVector item_ids      = items.column(0),
                statement_ids = items.column(1) - 1,
                dim_ids       = items.column(2),
                unique_items  = sort_unique(item_ids);
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

  return List::create(
    _["item_ids"]      = item_ids,
    _["statement_ids"] = statement_ids,
    _["dim_ids"]       = dim_ids,
    _["unique_items"]  = unique_items,
    _["n_items"]       = n_items,
    _["n_persons"]     = n_persons,
    _["n_dims"]        = n_dims_all,
    _["n_rows_items"]  = n_rows_items,
    _["n_cols_items"]  = n_cols_items
  );
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
  //  - picked_orders: matrix of picked_order_id for [people x items]

  // Value:
  //  - matrix of loglikelihoods for mupp rank stuff ... (potentially aggregated)

  // initialize parameters and return stuff
  List all_params  = initialize_mupp_l(thetas, params, items, picked_orders);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"],
      n_items                     = all_params["n_items"],
      n_rows_items                = all_params["n_rows_items"];
  IntegerVector unique_items      = all_params["unique_items"],
                item_ids          = all_params["item_ids"],
                statement_ids     = all_params["statement_ids"],
                dim_ids           = all_params["dim_ids"];

  // indicate return matrix and temporary vectors
  NumericMatrix loglik(Rf_allocMatrix(REALSXP, n_persons, n_items));
  NumericVector p         = no_init(n_persons);
  LogicalVector item_flag = no_init(n_rows_items);
  IntegerVector item_resp = no_init(n_persons);

  for(int item = 0; item < n_items; item++){

    // pulling out statements and dimensions for this particular item
    item_flag = (item_ids == unique_items[item]);
    item_resp = picked_orders.column(item);

    // calculating probability (all thetas, relevant params/dims, item orders)
    p         = p_mupp_rank_impl(thetas,
                                 select_rows(params, statement_ids[item_flag]),
                                 dim_ids[item_flag],
                                 item_resp).column(0);

    // calculating likelihood and adding elements to output matrix (fixing NAs)
    loglik.column(item) = ifelse(is_na(item_resp), NA_REAL, log(p));
  }

  return loglik;
}

//[[Rcpp::export]]
NumericMatrix lder1_mupp_rank_impl(const NumericMatrix & thetas,
                                   const NumericMatrix & params,
                                   const IntegerMatrix & items,
                                   const IntegerMatrix & picked_orders){

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for all params, in order)
  //  - items: matrix of [item, statement, dim] for all items, where
  //           statement is integer indicating the statement number and aligns
  //           with params
  //  - picked_orders: matrix of picked_order_id for [people x items]

  // Value:
  //  - matrix of loglikelihoods for mupp rank stuff ... (potentially aggregated)

  // initialize parameters and return stuff
  List all_params  = initialize_mupp_l(thetas, params, items, picked_orders);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"],
      n_items                     = all_params["n_items"],
      n_dims                      = all_params["n_dims"],
      n_rows_items                = all_params["n_rows_items"];
  IntegerVector unique_items      = all_params["unique_items"],
                item_ids          = all_params["item_ids"],
                statement_ids     = all_params["statement_ids"],
                dim_ids           = all_params["dim_ids"];

  // indicate return matrix and temporary vectors
  NumericVector p         = no_init(n_persons);
  NumericMatrix pder1(Rf_allocMatrix(REALSXP, n_persons, n_dims)),
                loglik(n_persons, n_dims);
  LogicalVector item_flag = no_init(n_rows_items);
  int dim;

  for(int item = 0; item < n_items; item++){

    // pulling out statements and dimensions for this particular item
    item_flag = (item_ids == unique_items[item]);

    // extracting params/dimensions/picked_orders for each item
    NumericMatrix item_par  = select_rows(params, statement_ids[item_flag]);
    IntegerVector item_dims = dim_ids[item_flag],
                  prev_dims(item_dims.size(), -1),
                  item_resp = picked_orders.column(item);

    // calculating probability/probability derivative for the chosen thing
    p     = p_mupp_rank_impl(thetas, item_par, item_dims, item_resp).column(0);
    pder1 = pder1_mupp_rank_impl(thetas, item_par, item_dims, item_resp)[0];

    // updating pder1 so that columns are divided by p
    for(int dim_id = 0; dim_id < item_dims.size(); dim_id++){
      dim = item_dims[dim_id] - 1;

      // make sure to only add column once :)
      if(any(prev_dims.begin(), prev_dims.end(), dim)){
        continue;
      } else{
        prev_dims[dim_id] = dim;
      }


      // adding elements to output matrix (ignoring elements if NA)
      for(int person = 0; person < n_persons; person++){
        if(item_resp[person] != NA_INTEGER){
          loglik[person + dim * n_persons] += pder1[person + dim * n_persons] / p[person];
        }
      }
    }
  }

  return loglik;
}

//[[Rcpp::export]]
NumericMatrix lder2_mupp_rank_impl(const NumericMatrix & thetas,
                                   const NumericMatrix & params,
                                   const IntegerMatrix & items,
                                   const IntegerMatrix & picked_orders){

  // Arguments:
  //  - thetas: matrix of persons x dims (for all dims)
  //  - params: matrix of dims x params (for all params, in order)
  //  - items: matrix of [item, statement, dim] for all items, where
  //           statement is integer indicating the statement number and aligns
  //           with params
  //  - picked_orders: matrix of picked_order_id for [people x items]

  // Value:
  //  - matrix of loglikelihoods for mupp rank stuff ... (potentially aggregated)

  // initialize parameters and return stuff
  List all_params  = initialize_mupp_l(thetas, params, items, picked_orders);

  // pull out useful parameter parts
  int n_persons                   = all_params["n_persons"],
      n_items                     = all_params["n_items"],
      n_dims                      = all_params["n_dims"],
      n_rows_items                = all_params["n_rows_items"];
  IntegerVector unique_items      = all_params["unique_items"],
                item_ids          = all_params["item_ids"],
                statement_ids     = all_params["statement_ids"],
                dim_ids           = all_params["dim_ids"];

  // number of pairs (off-diagonal)
  int n_pairs = n_dims * (n_dims - 1) / 2;

  // indicate return matrix and temporary vectors
  NumericVector p         = no_init(n_persons);
  NumericMatrix pder1(Rf_allocMatrix(REALSXP, n_persons, n_dims)),
                pder2(Rf_allocMatrix(REALSXP, n_persons, n_dims + n_pairs)),
                loglik(n_persons, n_dims + n_pairs);
  LogicalVector item_flag = no_init(n_rows_items);
  int dim1, dim2, dim12;
  double p1, dp1, dp2, dp12;

  for(int item = 0; item < n_items; item++){

    // pulling out statements and dimensions for this particular item
    item_flag = (item_ids == unique_items[item]);

    // extracting params/dimensions/picked_orders for each item
    NumericMatrix item_par  = select_rows(params, statement_ids[item_flag]);
    IntegerVector item_dims = dim_ids[item_flag],
                  prev_dims(n_dims + n_pairs, -1),
                  item_resp = picked_orders.column(item);

    // calculating probability/probability derivative/second derivatives for the chosen thing
    p     = p_mupp_rank_impl(thetas, item_par, item_dims, item_resp).column(0);
    pder1 = pder1_mupp_rank_impl(thetas, item_par, item_dims, item_resp)[0];
    pder2 = pder2_mupp_rank_impl(thetas, item_par, item_dims, item_resp)[0];

    // creating second derivative by combining p, pd1, and pd2
    for(int dim_id1 = 0; dim_id1 < item_dims.size(); dim_id1++){
      for(int dim_id2 = dim_id1; dim_id2 < item_dims.size(); dim_id2++){
        dim1  = item_dims[dim_id1] - 1;
        dim2  = item_dims[dim_id2] - 1;
        dim12 = find_crossprod_column(dim1, dim2, n_dims);

        // make sure to only add column once :)
        if(any(prev_dims.begin(), prev_dims.end(), dim12)){
          continue;
        } else{
          prev_dims[dim12] = dim12;
        }

        // determining/adding elements to output matrix
        for(int person = 0; person < n_persons; person++){

          // adding elements to output matrix (ignoring elements if NA)
          if(item_resp[person] != NA_INTEGER){

            // pulling out p, dp1, dp2, dp12 for current combination
            // (dp12 is cross-derivative/second derivative)
            p1   = p[person];
            dp1  = pder1[person + dim1  * n_persons];
            dp2  = pder1[person + dim2  * n_persons];
            dp12 = pder2[person + dim12 * n_persons];

            // adding elements to output matrix (ignoring elements if NA)
            loglik[person + dim12 * n_persons] += (p1 * dp12 - dp1 * dp2) / pow(p1, 2);
          }
        }
      }
    }
  }

  return loglik;
}
