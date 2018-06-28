#ifndef GGUM_PROBS_H
#define GGUM_PROBS_H

#include <Rcpp.h>
double exp_ggum(const double theta,
                const double alpha,
                const double delta,
                double tau,
                int    exp_mult = 1);

Rcpp::NumericVector q_ggum(const SEXP & thetas,
                           const SEXP & params);

Rcpp::NumericMatrix q_ggum_all(const Rcpp::NumericMatrix & thetas,
                               const Rcpp::NumericMatrix & params);

Rcpp::NumericMatrix p_ggum_all(const Rcpp::NumericMatrix & thetas,
                               const Rcpp::NumericMatrix & params);

Rcpp::NumericVector pder1_ggum(const SEXP & thetas,
                               const SEXP & params);

Rcpp::NumericMatrix pder1_ggum_all(const Rcpp::NumericMatrix & thetas,
                                   const Rcpp::NumericMatrix & params);
#endif
