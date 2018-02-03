/* These functions are designed to calculate (mupp)
 * - probabilities
 * - first derviatives of probabilties
 * - second derivatives of probabilities
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector p_ggum(NumericVector thetas, NumericVector params){

  // Arguments:
  //  - thetas: a vector of thetas across all people
  //  - params: a vector of params for one item [alpha, delta, tau]
  // Value:
  //  - P(Z_s = 1 | thetas)

  // declare parameters
  double alpha = params[0];
  double delta = params[1];
  double tau   = params[2];

  // calculate exponent expressions
  NumericVector diff   = thetas - delta;
  NumericVector exp_12 = exp(alpha * (1 * diff - tau)) +
                         exp(alpha * (2 * diff - tau));
  NumericVector exp_3  = exp(alpha * (3 * diff));

  // calculate probabilities
  NumericVector probs  = exp_12 / (1 + exp_12 + exp_3);

  return probs;

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

    // set probs to item column of matrix
    probs(_, item) = (p_s1 * p_t0) / (p_s1 * p_t0 + (1 - p_s1) * (1 - p_t0));
  }

  return probs;

}


// R Test Code
/*** R
thetas   <- seq(-3, 3, length.out = 100000)
params_1 <- c(1.1, -2.3, -3.1)
params_2 <- c(1.4, 1.1, -2.4)

library(microbenchmark)

microbenchmark({
p_s1 <- p_ggum(thetas, params_1)
p_t1 <- p_ggum(thetas, params_2)
p_s0 <- 1 - p_s1
p_t0 <- 1 - p_t1

out1  <- p_s1 * p_t0 / ((p_s1 * p_t0) + (p_s0 * p_t1))
},
{
out2  <- p_mupp(thetas, thetas, rbind(params_1), rbind(params_2))
})
*/
