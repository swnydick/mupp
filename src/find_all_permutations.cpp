#include <Rcpp.h>
using namespace Rcpp;

//' Find All Permutations of Consecutive Integers
//'
//' Given an consecutive integer vector of length n, find all permutations of
//' that vector.
//'
//' @param n an integer greater than 0
//' @param init an integer indicating the initial starting value for the set
//'        of integers included in the permutation. See Details.
//'
//' @return A matrix of size n! x n, where each row is a unique permutation
//'
//' @details \code{init} is useful for indexing C vs R code. If \code{init = 0},
//'          then the indices will work with 0 indexed languages, such as C or
//'          Python. If \code{init = 1}, then the indices will work with 1 indexed
//'          languages, such as R.
//'
//' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
//'
//' @export
// [[Rcpp::export]]
IntegerMatrix find_all_permutations(int n,
                                    int init = 0){

  // Arguments:
  //  - n: an integer greater than or equal to 1
  // Value:
  //  - permutation matrix in reasonable order

  // Argument Checks
  if(n < 1){
    stop("n must be a scalar greater than 1");
  }

  // declaring helper vectors
  IntegerVector n_r(1, n);

  // declaring parameter parts
  int n_rows   = factorial(n_r)[0],
      n_cols   = n,
      i        = 0;

  // indicating return matrix and temporary storage vector
  IntegerVector cur_perm = seq_len(n_cols) - (1 - init);
  IntegerMatrix all_perms(n_rows, n_cols);

  // add permutation to all_perms vector, then
  // use iterators to find next permutation
  do{
    all_perms(i, _) = cur_perm;
    i              += 1;
  } while (std::next_permutation(cur_perm.begin(), cur_perm.end()));

  return all_perms;
}


//' Find Column for Cross-Product
//'
//' Given two dimension and the total number of dimensions, find the column
//' indexing the appropriate cross-product.
//'
//' @param dim1 the first dimension
//' @param dim2 the second dimension
//' @param n_dims the total number of dimensions
//' @param init an integer indicating the initial starting value for the set
//'        of integers included in the permutation. See Details.
//'
//' @return The index of the cross-product column.
//'
//' @details \code{init} is useful for indexing C vs R code. If \code{init = 0},
//'          then the indices will work with 0 indexed languages, such as C or
//'          Python. If \code{init = 1}, then the indices will work with 1 indexed
//'          languages, such as R.
//'
//'          This function assumes that the first n_dims columns are the dims,
//'          the next n_dims - 1 columns are the cross-products of dimension 1 with
//'          the remaining dims, the next n_dims - 2 columns are the cross-products
//'          of dimension 2 with the remaining dims > 2, etc.
//'
//' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
//'
//' @export
// [[Rcpp::export]]
int find_crossprod_column(int dim1,
                          int dim2,
                          int n_dims,
                          int init = 0) {

  // integer to store cross-product column
  int dim12 = n_dims - 1;
  dim1 -= init;
  dim2 -= init;

  // first  n_dim columns are dim1, dim2, ...
  // second n_dim - 1 columns are dim12, dim13, dim14, ...
  // third  n_dim - 1 columns are dim23, dim24, ...

  // argument checks //
  if(n_dims <= std::max(dim1, dim2)){
    stop("n_dim is not at least the number of dimensions in dim1_col OR dim2_col");
  }

  // fixing arguments:
  //  - if dim1 > dim2, swap!
  //  - if dim1 or dim2 < 0, bad!
  //  - if dim1 = dim2, return -1 (no cross-product)
  if(dim1 > dim2){
    int dim3 = dim1;
    dim1 = dim2;
    dim2 = dim3;
  } else if(dim1 < 0 | dim2 < 0){
    stop("dim1 and dim2 must greater than init");
  } else if(dim1 == dim2){
    dim12 = -1;
  } else{

    // from 0 to dim1, add the next possible unique combs (with dim)
    for(int dim = 0; dim < dim1; dim++){
      dim12 += n_dims - dim - 1;
    }

    // add the remainder
    dim12 += dim2 - dim1;
  }

  return dim12 + init;
}
