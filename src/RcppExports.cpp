// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// cp_gibbsC
NumericMatrix cp_gibbsC(int length_of_chains, int from_point, int a, int b, int x_range, float mu_x, float mu_y);
RcppExport SEXP _StatComp21047_cp_gibbsC(SEXP length_of_chainsSEXP, SEXP from_pointSEXP, SEXP aSEXP, SEXP bSEXP, SEXP x_rangeSEXP, SEXP mu_xSEXP, SEXP mu_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type length_of_chains(length_of_chainsSEXP);
    Rcpp::traits::input_parameter< int >::type from_point(from_pointSEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type x_range(x_rangeSEXP);
    Rcpp::traits::input_parameter< float >::type mu_x(mu_xSEXP);
    Rcpp::traits::input_parameter< float >::type mu_y(mu_ySEXP);
    rcpp_result_gen = Rcpp::wrap(cp_gibbsC(length_of_chains, from_point, a, b, x_range, mu_x, mu_y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_StatComp21047_cp_gibbsC", (DL_FUNC) &_StatComp21047_cp_gibbsC, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_StatComp21047(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
