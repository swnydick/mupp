#ifndef MUPP_UTILITIES_H
#define MUPP_UTILITIES_H

#include <Rcpp.h>
Rcpp::IntegerMatrix find_all_permutations(int n, int init = 0);

int find_crossprod_column(int dim1, int dim2, int n_dims, int init = 0);

#endif
