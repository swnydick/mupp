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
  int n_rows   = factorial(n_r)[0];
  int n_cols   = n;
  int i        = 0;

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

