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
 *   - find_all_permutations (permutation matrix)
 */
NumericMatrix select_cols(NumericMatrix X,
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

/* GGUM STUFF
 *  - exp_ggum
 *  - p_ggum/q_ggum_all
 *  - p_der1_theta_ggum
 */
NumericVector exp_ggum(NumericVector thetas,
                       NumericVector params,
                       int exp_mult = 1) {

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

NumericVector q_ggum(NumericVector thetas,
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

NumericMatrix q_ggum_all(NumericMatrix thetas,
                         NumericMatrix params) {

  // Arguments:
  //  - thetas: a matrix of thetas across all people/dims
  //  - params: a matrix of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z = 0 | thetas) across all dimensions

  // declare number of persons/dimensions
  int n_persons = thetas.nrow();
  int n_dims    = params.nrow();

  // indicating return vector
  NumericMatrix probs(n_persons, n_dims);

  // calculating probabilities across all dimensions
  for(int dim = 0; dim < n_dims; dim++){
    probs(_, dim) = q_ggum(thetas(_, dim), params(dim, _));
  }

  return probs;
}

NumericVector pder1_theta_ggum(NumericVector thetas,
                               NumericVector params) {

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
NumericMatrix pder1_theta_ggum_all(NumericMatrix thetas,
                                   NumericMatrix params) {

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

NumericVector p_mupp_pick0(NumericMatrix Q,
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

NumericVector p_mupp_pick1(NumericMatrix Q,
                           int picked_dim = 1) {

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

// [[Rcpp::export]]
NumericMatrix p_mupp_rank_impl(NumericMatrix thetas,
                               NumericMatrix params,
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

  // fixing dims if it is NA or missing
  if(any(is_na(dims))){
    dims = seq_len(n_params);
  }

  if(any(is_na(picked_order_id))){
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
  IntegerMatrix picked_orders = find_all_permutations(n_dims_item);
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
    NumericVector Q_person(n_params);
    IntegerVector picked_order(n_params);

    // for each person, rearrange Q so that [1, 2, 3, ...] is PICKED order ...
    for(int person = 0; person < n_persons; person++){
      picked_order = picked_orders(picked_order_id[person], _);
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


/*** R
rm(list = ls())

# specifying parameters
params_1 <- c(2, -2, -1)
params_2 <- c(2, 0, 0)
params_3 <- c(2, 1, -1)
params_4 <- c(2, .1, .1)
params_5 <- c(2, -3, 3)
params_6 <- c(1, .8, 0)
params_7 <- c(1, -2, -.2)

# generating parameters/thetas based on specification
params   <- do.call(what = rbind,
                    args = lapply(ls(pattern = "params\\_",
                                     envir   = .GlobalEnv),
                                  FUN = get))
thetas   <- seq(-3, 3, length.out = 1000)
thetas   <- do.call(what = cbind,
                    args = lapply(1:nrow(params),
                                  FUN = function(x) thetas))

# calculating
out      <- mupp:::p_mupp_rank_impl(thetas, params)

# restructuring
comb     <- apply(mupp:::find_all_permutations(nrow(params), 1),
                  MARGIN   = 1,
                  FUN      = paste,
                  collapse = "-")
out      <- setNames(object = as.data.frame(out),
                     nm     = comb)
out      <- cbind(theta = thetas[ , 1], out)
out      <- reshape2::melt(out,
                           id.vars       = "theta",
                           variable.name = "combination",
                           value.name    = "probability")

# plotting
library(ggplot2)
g <- ggplot(out, aes(x        = theta,
                     y        = probability,
                     color    = combination)) +
     geom_line(size = 1) +
     theme_minimal() +
     guides(color = FALSE)
print(g)

# derivatives??
library(numDeriv)

p_mupp_rank <- function(thetas, params, rank_index = 1){

  mupp:::p_mupp_rank_impl(rbind(thetas), rbind(params))[ , rank_index]

}
*/
