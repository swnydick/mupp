#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix select_cols(const NumericMatrix & X,
                          const IntegerVector & ind) {

  // Arguments:
  //  - X:   a matrix to select columns
  //  - ind: an integer vector of columns to select

  // capturing constants
  int out_cols = ind.size(),
      out_rows = X.nrow(),
      inp_cols = X.ncol();
  IntegerVector inp_compare = seq_len(inp_cols) - 1;

  // - if all ind is equal to the number of columns of X, return X
  // - otherwise, re-arrange
  if((out_cols == inp_compare.size()) & is_true(all(ind == inp_compare))){
    return clone(X);
  } else{
    NumericMatrix Y(Rf_allocMatrix(REALSXP, out_rows, out_cols));

    for(int j = 0; j < out_cols; j++){
      for(int i = 0; i < out_rows; i++){
        Y[i + j * out_rows] = X[i + ind[j] * out_rows];
      }
    }

    return Y;
  }
}

NumericMatrix select_rows(const NumericMatrix & X,
                          const IntegerVector & ind) {

  // Arguments:
  //  - X:   a matrix to select rows
  //  - ind: an integer vector of rows to select

  // capturing constants
  int out_rows = ind.size(),
      out_cols = X.ncol(),
      inp_rows = X.nrow();
  IntegerVector inp_compare = seq_len(inp_rows) - 1;

  // - if all ind is equal to the number of columns of X, return X
  // - otherwise, re-arrange
  if((out_rows == inp_compare.size()) & is_true(all(ind == inp_compare))){
    return clone(X);
  } else{
    NumericMatrix Y(Rf_allocMatrix(REALSXP, out_rows, out_cols));

    for(int j = 0; j < out_cols; j++){
      for(int i = 0; i < out_rows; i++){
        Y[i + j * out_rows] = X[ind[i] + j * inp_rows];
      }
    }

    return Y;
  }
}


// arrange matrix so that every row is sorted by picked order
void arrange_by_picked(NumericMatrix M,
                       const IntegerVector & indices,
                       const IntegerMatrix & orders,
                       bool reverse = false) {

  // declare number of rows/columns
  int n_dims    = M.ncol(),
      n_persons = M.nrow(),
      n_orders  = orders.ncol();

  // temporary vectors to store the matrices and picked order
  NumericVector M_person   = no_init(n_dims);
  IntegerMatrix all_orders = transpose(orders);

  // indices to reverse the ordering (by row) if required)
  if(reverse){
    for(int col = 0; col < n_orders; col++){
      IntegerVector col_order   = all_orders.column(col);
      IntegerVector match_order = clone(col_order).sort();
      all_orders.column(col)    = match(match_order, col_order) - 1;
    }
  }

  // for each person, rearrange M so that [1, 2, 3] is the PICKED order ...
  for(int person = 0; person < n_persons; person++){

    // pull out appropriate row of matrix (required due to swapping)
    M_person = M.row(person);

    // put the ordered element in the appropriate place in the matrix (ignoring NAs)
    for(int i = 0; i < n_orders; i++){
      if(indices[person] != NA_INTEGER){
        M[person + n_persons * i] = M_person[all_orders[i + n_orders * indices[person]]];
      }
    }
  }
}
