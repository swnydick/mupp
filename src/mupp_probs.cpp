/* These functions are designed to calculate (mupp)
 * - probabilities
 * - first derviatives of probabilties
 * - second derivatives of probabilities
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
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
NumericVector p_ggum(NumericVector thetas, NumericVector params){

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
NumericVector pder1_theta_ggum(NumericVector thetas, NumericVector params){

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
NumericMatrix p_mupp(NumericVector thetas_s, NumericVector thetas_t,
                     NumericMatrix params_s, NumericMatrix params_t) {

  // Arguments:
  //  - thetas_s/thetas_t: vectors of thetas across all people
  //  - params_s/params_t: matrices of params for all items
  // Value:
  //  - P(s > t)(theta_ss, thetas_t) for all items

  // declare number of items and persons
  int n_items   = std::max(params_s.nrow(), params_t.nrow());
  int n_persons = std::max(thetas_s.size(), thetas_t.size());
  int n_params  = std::max(params_s.ncol(), params_t.ncol());

  // indicate temporary storage vectors
  NumericVector params_si(n_params),
                params_ti(n_params),
                p_s1(n_persons),
                p_t0(n_persons);

  // indicate return matrix
  NumericMatrix probs(n_persons, n_items);

  // cycle through and add vectors to return matrix
  for(int item = 0; item < n_items; item++){

    // params for this specific item
    params_si = params_s(item, _);
    params_ti = params_t(item, _);

    // probability of s/t given ggum model (can use other models?)
    p_s1 = p_ggum(thetas_s, params_si);
    p_t0 = 1 - p_ggum(thetas_t, params_ti);

    // set probs to item column of probability matrix
    probs(_, item) = (p_s1 * p_t0) / (p_s1 * p_t0 + (1 - p_s1) * (1 - p_t0));
  }

  return probs;

}

// [[Rcpp::export]]
NumericMatrix pder1_thetas_mupp_(NumericVector thetas_s, NumericVector thetas_t,
                                 NumericMatrix params_s, NumericMatrix params_t,
                                 bool type_s = true) {

  // Arguments:
  //  - thetas_s/thetas_t: vectors of thetas across all people
  //  - params_s/params_t: matrices of params for all items
  // Value:
  //  - dp(s > t)/dtheta(theta_ss, thetas_t) for all items

  // declare number of items and persons
  int n_items   = std::max(params_s.nrow(), params_t.nrow());
  int n_persons = std::max(thetas_s.size(), thetas_t.size());
  int n_params  = std::max(params_s.ncol(), params_t.ncol());

  // indicate temporary storage vectors
  NumericVector params_si(n_params),
                params_ti(n_params),
                p_s1(n_persons),
                p_t0(n_persons),
                dp_s1(n_persons),
                dp_t0(n_persons),
                denom(n_persons);

  // indicate return matrix
  NumericMatrix dprobs(n_persons, n_items);

  // cycle through and add vectors to return matrix
  for(int item = 0; item < n_items; item++){

    // params for this specific item
    params_si = params_s(item, _);
    params_ti = params_t(item, _);

    // probability of s/t given ggum model (can use other models?)
    p_s1  = p_ggum(thetas_s, params_si);       // A
    p_t0  = 1 - p_ggum(thetas_t, params_ti);   // B

    // derivatives of s/t given ggum model (can use other models?)
    dp_s1 =  pder1_theta_ggum(thetas_s, params_si);  // A'
    dp_t0 = -pder1_theta_ggum(thetas_s, params_ti);  // B'

    // denominator of the model
    denom = p_s1 * p_t0 + (1 - p_s1) * (1 - p_t0);

    // set dprobs to item column of derivative matrix
    if(type_s){
      dprobs(_, item) = ((dp_s1 * p_t0) * (denom + p_s1)) / pow(denom, 2);
    } else{
      dprobs(_, item) = ((dp_t0 * p_s1) * (denom + p_t0)) / pow(denom, 2);
    }
  }

  return dprobs;

}


// R Test Code
/*** R
thetas   <- seq(-3, 3, length.out = 100000)
params_1 <- c(1.1, -2.3, -3.1)
params_2 <- c(1.4, 1.1, -2.4)

library(microbenchmark)

microbenchmark(v1 = {
p_s1 <- p_ggum(thetas, params_1)
p_t1 <- p_ggum(thetas, params_2)
p_s0 <- 1 - p_s1
p_t0 <- 1 - p_t1

out1  <- p_s1 * p_t0 / ((p_s1 * p_t0) + (p_s0 * p_t1))
},
v2 = {
out2  <- p_mupp(thetas, thetas, rbind(params_1), rbind(params_2))
})
*/
