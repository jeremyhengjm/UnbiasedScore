// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// multinomial_resampling
IntegerVector multinomial_resampling(const NumericVector& weights, int ndraws, const NumericVector& rand);
RcppExport SEXP _UnbiasedScore_multinomial_resampling(SEXP weightsSEXP, SEXP ndrawsSEXP, SEXP randSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type rand(randSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomial_resampling(weights, ndraws, rand));
    return rcpp_result_gen;
END_RCPP
}
// maximal_rejection_sampling
IntegerVector maximal_rejection_sampling(const IntegerVector& ancestors1, const NumericVector& weights1, const NumericVector& weights2);
RcppExport SEXP _UnbiasedScore_maximal_rejection_sampling(SEXP ancestors1SEXP, SEXP weights1SEXP, SEXP weights2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type ancestors1(ancestors1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights1(weights1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights2(weights2SEXP);
    rcpp_result_gen = Rcpp::wrap(maximal_rejection_sampling(ancestors1, weights1, weights2));
    return rcpp_result_gen;
END_RCPP
}
// maximal_maximal_multinomial_resampling
IntegerMatrix maximal_maximal_multinomial_resampling(const NumericVector& weights1, const NumericVector& weights2, const NumericVector& weights3, const NumericVector& weights4, int ndraws, int residualtype);
RcppExport SEXP _UnbiasedScore_maximal_maximal_multinomial_resampling(SEXP weights1SEXP, SEXP weights2SEXP, SEXP weights3SEXP, SEXP weights4SEXP, SEXP ndrawsSEXP, SEXP residualtypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type weights1(weights1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights2(weights2SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights3(weights3SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights4(weights4SEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    Rcpp::traits::input_parameter< int >::type residualtype(residualtypeSEXP);
    rcpp_result_gen = Rcpp::wrap(maximal_maximal_multinomial_resampling(weights1, weights2, weights3, weights4, ndraws, residualtype));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_module_tree();

static const R_CallMethodDef CallEntries[] = {
    {"_UnbiasedScore_multinomial_resampling", (DL_FUNC) &_UnbiasedScore_multinomial_resampling, 3},
    {"_UnbiasedScore_maximal_rejection_sampling", (DL_FUNC) &_UnbiasedScore_maximal_rejection_sampling, 3},
    {"_UnbiasedScore_maximal_maximal_multinomial_resampling", (DL_FUNC) &_UnbiasedScore_maximal_maximal_multinomial_resampling, 6},
    {"_rcpp_module_boot_module_tree", (DL_FUNC) &_rcpp_module_boot_module_tree, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_UnbiasedScore(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
