// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bspline_basis
double bspline_basis(int j, int m, NumericVector knots, double x);
RcppExport SEXP _smooth_commutability_bspline_basis(SEXP jSEXP, SEXP mSEXP, SEXP knotsSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(bspline_basis(j, m, knots, x));
    return rcpp_result_gen;
END_RCPP
}
// bspline_basis_matrix
NumericMatrix bspline_basis_matrix(NumericVector x, NumericVector knots, int m);
RcppExport SEXP _smooth_commutability_bspline_basis_matrix(SEXP xSEXP, SEXP knotsSEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type knots(knotsSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(bspline_basis_matrix(x, knots, m));
    return rcpp_result_gen;
END_RCPP
}
// df_2
double df_2(NumericVector df_grid, NumericVector cv, double df_0, double D, double tol, bool silence);
RcppExport SEXP _smooth_commutability_df_2(SEXP df_gridSEXP, SEXP cvSEXP, SEXP df_0SEXP, SEXP DSEXP, SEXP tolSEXP, SEXP silenceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type df_grid(df_gridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cv(cvSEXP);
    Rcpp::traits::input_parameter< double >::type df_0(df_0SEXP);
    Rcpp::traits::input_parameter< double >::type D(DSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< bool >::type silence(silenceSEXP);
    rcpp_result_gen = Rcpp::wrap(df_2(df_grid, cv, df_0, D, tol, silence));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _smooth_commutability_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// simulate_eqa_data_custom_cpp
List simulate_eqa_data_custom_cpp(List parameters, int type, bool AR);
RcppExport SEXP _smooth_commutability_simulate_eqa_data_custom_cpp(SEXP parametersSEXP, SEXP typeSEXP, SEXP ARSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type parameters(parametersSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< bool >::type AR(ARSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_eqa_data_custom_cpp(parameters, type, AR));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_smooth_commutability_bspline_basis", (DL_FUNC) &_smooth_commutability_bspline_basis, 4},
    {"_smooth_commutability_bspline_basis_matrix", (DL_FUNC) &_smooth_commutability_bspline_basis_matrix, 3},
    {"_smooth_commutability_df_2", (DL_FUNC) &_smooth_commutability_df_2, 6},
    {"_smooth_commutability_rcpp_hello_world", (DL_FUNC) &_smooth_commutability_rcpp_hello_world, 0},
    {"_smooth_commutability_simulate_eqa_data_custom_cpp", (DL_FUNC) &_smooth_commutability_simulate_eqa_data_custom_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_smooth_commutability(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
