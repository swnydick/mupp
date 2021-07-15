#include <Rcpp.h>
using namespace Rcpp;

// determine rank of vector (in utilities???)
IntegerVector rank_ints(IntegerVector x){

  // clone the vector and sort it
  IntegerVector sorted = clone(x).sort();

  // determine which variables are sorted
  return match(x, sorted);
}

//' Find All Permutations of Consecutive Integers
//'
//' Given an consecutive integer vector of length n, find all permutations of
//' that vector.
//'
//' @param n an integer greater than 0
//' @param init an integer indicating the initial starting value for the set
//'        of integers included in the permutation. See Details.
//' @param index an integer between 1 and n indicatng which permutation to
//'        select (using R indexing rather than C indexing)
//' @param order an integer vector indicating the response/permutation order
//'        that should be checked against
//'
//' @return For `find_all_permutations`, a A matrix of size n! x n, where each
//'         row is a unique permutation. For `find_permutation_order`, the
//'         index row into the permutation matrix (where "index" uses R indexing).
//'         For `find_permutation_index` the index of the permutation matrix
//'         that leads to a particular order.
//'
//' @note `find_permutation_order` and `find_permutation_index` both use R indexing
//'       for determining the row in the permutation matrix and are almost
//'       inverses of each other. The reason they are "almost" inverses and not
//'       true inverses is because `find_permutation_index` needs only an
//'       integer vector and will standardize that vector to be between 1 and
//'       n. Therefore, the order is unique only up to a monotonic transformation.
//'
//' @details \code{init} is useful for indexing C vs R code. If \code{init = 0},
//'          then the indices will work with 0 indexed languages, such as C or
//'          Python. If \code{init = 1}, then the indices will work with 1 indexed
//'          languages, such as R.
//'
//' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
//' @name permutation

//' @rdname permutation
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

// saved permutations (so we don't have to recreate this all of the time)
List saved_permutations = List::create(find_all_permutations(1),
                                       find_all_permutations(2),
                                       find_all_permutations(3),
                                       find_all_permutations(4),
                                       find_all_permutations(5),
                                       find_all_permutations(6),
                                       find_all_permutations(7));

// extract save permutation (if we have it stored)
IntegerMatrix extract_permutations(int n,
                                   int init = 0){

  if((n <= saved_permutations.size()) & (init == 0)){
    return saved_permutations[n - 1];
  } else{
    return find_all_permutations(n, init);
  }
}

//' @rdname permutation
//' @export
// [[Rcpp::export]]
IntegerVector find_permutation_order(int n,
                                     int index = 1,
                                     int init  = 0){

  // pull out the appropriate permutation set
  IntegerMatrix picked_orders = extract_permutations(n, init);

  // Argument Checks
  if(index < 1 | index > picked_orders.nrow()){
    stop("index must be between 1 and n!");
  }

  // pull out the appropriate permutation
  return picked_orders(index - 1, _);
}

//' @rdname permutation
//' @export
// [[Rcpp::export]]
IntegerVector find_permutation_index(IntegerVector order){

  // pull out the size of the vector
  int n = order.length();

  // determine whether the vector has any problems (NA or duplicated values)
  bool has_problem = is_true(any(is_na(order) | duplicated(order)));

  // end early and return something sensible if short length or NAs
  if(n == 1){
    return {1};
  } else if((n == 0) | has_problem){
    return {NA_INTEGER};
  }

  // pull out the appropriate permutation set and reorder vector
  IntegerMatrix picked_orders = extract_permutations(n);
  order = rank_ints(order) - 1;

  // iteratively determine if any of the orders are the same
  for(int i = 0; i < picked_orders.nrow(); i++){
    if(is_true(all(order == picked_orders(i, _)))){
      return {i + 1};
    }
  }

  // return NA if something screwed up somewhere above
  return {NA_INTEGER};
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
//' @seealso \code{\link{find_crossprod_dims}}
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
    stop("n_dim is not at least the number of dimensions in dim1 OR dim2");
  }

  // fixing arguments:
  //  - if dim1 > dim2, swap!
  //  - if dim1 or dim2 < 0, bad!
  //  - if dim1 = dim2, return the element of the diagonal
  if((dim1 < 0) | (dim2 < 0)){
    stop("dim1 and dim2 must greater than init");
  } else if(dim1 > dim2){
    int dim3 = dim1;
    dim1 = dim2;
    dim2 = dim3;
  } else if(dim1 == dim2){
    dim12 = dim1;
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


//' Find Dimensions given Cross-Product Column
//'
//' Given the column indexing the appropriate cross-product and the total number
//' of dimensions, find the two dimensions leading to the cross-product.
//'
//' @param dim12 the column of the matrix or data.frame
//' @param n_dims the total number of dimensions
//' @param init an integer indicating the initial starting value for the set
//'        of integers included in the permutation. See Details.
//'
//' @return A vector/array indicating the cross-product dimensions.
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
//' @seealso \code{\link{find_crossprod_column}}
//'
//' @author Steven Nydick, \email{steven.nydick@@kornferry.com}
//'
//' @export
// [[Rcpp::export]]
NumericVector find_crossprod_dims(int dim12,
                                  int n_dims,
                                  int init = 0) {

  // integer to store dimension columns
  int dim1, dim2;
  dim12 -= init;

  // first  n_dim columns are dim1, dim2, ...
  // second n_dim - 1 columns are dim12, dim13, dim14, ...
  // third  n_dim - 1 columns are dim23, dim24, ...

  // fixing arguments:
  //  - if dim12 < 0, bad!
  //  - if dim12 < n_dims, return the element of the diagonal
  if(dim12 < 0){
    stop("dim12 must greater than init");
  } else if(dim12 < n_dims){
    dim1 = dim12;
    dim2 = dim12;
  } else{

    // first, subtract the diagonal
    dim12 -= n_dims;

    // then subtract the dimension from dim12 until we go below 0
    for(dim1 = 0; dim1 < n_dims; dim1++){
      dim12 -= n_dims - dim1 - 1;

      if(dim12 < 0){
        break;
      }
    }

    // add the remainder
    dim2 = n_dims + dim12;
  }

  // argument checks //
  if(n_dims <= dim2){
    stop("n_dim is not at least the number of dimensions in dim1_col OR dim2_col");
  }

  return NumericVector::create(dim1 + init,
                               dim2 + init);
}
