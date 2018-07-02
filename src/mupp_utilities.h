#ifndef MUPP_UTILITIES_H
#define MUPP_UTILITIES_H

#include <Rcpp.h>
Rcpp::IntegerMatrix find_all_permutations(int n, int init = 0);

int find_crossprod_column(int dim1, int dim2, int n_dims, int init = 0);
Rcpp::NumericVector find_crossprod_dims(int dim12, int n_dims, int init = 0);

Rcpp::NumericMatrix select_cols(const Rcpp::NumericMatrix & X,
                                const Rcpp::IntegerVector & ind);
Rcpp::NumericMatrix select_rows(const Rcpp::NumericMatrix & X,
                                const Rcpp::IntegerVector & ind);

void arrange_by_picked(Rcpp::NumericMatrix M,
                       const Rcpp::IntegerVector & indices,
                       const Rcpp::IntegerMatrix & orders,
                       bool reverse = false);

#endif
