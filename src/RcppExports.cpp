// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// start_profiler
SEXP start_profiler(SEXP str);
RcppExport SEXP _mupp_start_profiler(SEXP strSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type str(strSEXP);
    rcpp_result_gen = Rcpp::wrap(start_profiler(str));
    return rcpp_result_gen;
END_RCPP
}
// stop_profiler
SEXP stop_profiler();
RcppExport SEXP _mupp_stop_profiler() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(stop_profiler());
    return rcpp_result_gen;
END_RCPP
}
// find_all_permutations
IntegerMatrix find_all_permutations(int n, int init);
RcppExport SEXP _mupp_find_all_permutations(SEXP nSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(find_all_permutations(n, init));
    return rcpp_result_gen;
END_RCPP
}
// extract_permutations
IntegerMatrix extract_permutations(int n, int init);
RcppExport SEXP _mupp_extract_permutations(SEXP nSEXP, SEXP initSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type init(initSEXP);
    rcpp_result_gen = Rcpp::wrap(extract_permutations(n, init));
    return rcpp_result_gen;
END_RCPP
}
// pder1_mupp_rank1
NumericMatrix pder1_mupp_rank1(const NumericMatrix& P, const NumericMatrix& dP, const IntegerVector order);
RcppExport SEXP _mupp_pder1_mupp_rank1(SEXP PSEXP, SEXP dPSEXP, SEXP orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type P(PSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type dP(dPSEXP);
    Rcpp::traits::input_parameter< const IntegerVector >::type order(orderSEXP);
    rcpp_result_gen = Rcpp::wrap(pder1_mupp_rank1(P, dP, order));
    return rcpp_result_gen;
END_RCPP
}
// p_mupp_rank_impl
NumericMatrix p_mupp_rank_impl(const NumericMatrix& thetas, const NumericMatrix& params, IntegerVector dims, IntegerVector picked_order_id);
RcppExport SEXP _mupp_p_mupp_rank_impl(SEXP thetasSEXP, SEXP paramsSEXP, SEXP dimsSEXP, SEXP picked_order_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type picked_order_id(picked_order_idSEXP);
    rcpp_result_gen = Rcpp::wrap(p_mupp_rank_impl(thetas, params, dims, picked_order_id));
    return rcpp_result_gen;
END_RCPP
}
// pder1_mupp_rank_impl
List pder1_mupp_rank_impl(NumericMatrix& thetas, const NumericMatrix& params, IntegerVector dims, IntegerVector picked_order_id);
RcppExport SEXP _mupp_pder1_mupp_rank_impl(SEXP thetasSEXP, SEXP paramsSEXP, SEXP dimsSEXP, SEXP picked_order_idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type dims(dimsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type picked_order_id(picked_order_idSEXP);
    rcpp_result_gen = Rcpp::wrap(pder1_mupp_rank_impl(thetas, params, dims, picked_order_id));
    return rcpp_result_gen;
END_RCPP
}
// loglik_mupp_rank_impl
NumericMatrix loglik_mupp_rank_impl(const NumericMatrix& thetas, const NumericMatrix& params, const IntegerMatrix& items, const IntegerMatrix& picked_orders);
RcppExport SEXP _mupp_loglik_mupp_rank_impl(SEXP thetasSEXP, SEXP paramsSEXP, SEXP itemsSEXP, SEXP picked_ordersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type items(itemsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type picked_orders(picked_ordersSEXP);
    rcpp_result_gen = Rcpp::wrap(loglik_mupp_rank_impl(thetas, params, items, picked_orders));
    return rcpp_result_gen;
END_RCPP
}
// lder1_mupp_rank_impl
List lder1_mupp_rank_impl(NumericMatrix& thetas, const NumericMatrix& params, const IntegerMatrix& items, const IntegerMatrix& picked_orders);
RcppExport SEXP _mupp_lder1_mupp_rank_impl(SEXP thetasSEXP, SEXP paramsSEXP, SEXP itemsSEXP, SEXP picked_ordersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type thetas(thetasSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type items(itemsSEXP);
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type picked_orders(picked_ordersSEXP);
    rcpp_result_gen = Rcpp::wrap(lder1_mupp_rank_impl(thetas, params, items, picked_orders));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mupp_start_profiler", (DL_FUNC) &_mupp_start_profiler, 1},
    {"_mupp_stop_profiler", (DL_FUNC) &_mupp_stop_profiler, 0},
    {"_mupp_find_all_permutations", (DL_FUNC) &_mupp_find_all_permutations, 2},
    {"_mupp_extract_permutations", (DL_FUNC) &_mupp_extract_permutations, 2},
    {"_mupp_pder1_mupp_rank1", (DL_FUNC) &_mupp_pder1_mupp_rank1, 3},
    {"_mupp_p_mupp_rank_impl", (DL_FUNC) &_mupp_p_mupp_rank_impl, 4},
    {"_mupp_pder1_mupp_rank_impl", (DL_FUNC) &_mupp_pder1_mupp_rank_impl, 4},
    {"_mupp_loglik_mupp_rank_impl", (DL_FUNC) &_mupp_loglik_mupp_rank_impl, 4},
    {"_mupp_lder1_mupp_rank_impl", (DL_FUNC) &_mupp_lder1_mupp_rank_impl, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mupp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
