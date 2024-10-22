// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// smat_to_mat
Eigen::MatrixXd smat_to_mat(const Eigen::SparseMatrix<double>& sparseMat, int k, bool equal_spaced);
RcppExport SEXP _trendfilter_smat_to_mat(SEXP sparseMatSEXP, SEXP kSEXP, SEXP equal_spacedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::SparseMatrix<double>& >::type sparseMat(sparseMatSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type equal_spaced(equal_spacedSEXP);
    rcpp_result_gen = Rcpp::wrap(smat_to_mat(sparseMat, k, equal_spaced));
    return rcpp_result_gen;
END_RCPP
}
// configure_denseD_test
Rcpp::List configure_denseD_test(Rcpp::NumericVector x, int k);
RcppExport SEXP _trendfilter_configure_denseD_test(SEXP xSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(configure_denseD_test(x, k));
    return rcpp_result_gen;
END_RCPP
}
// linear_single_solve_test
Eigen::VectorXd linear_single_solve_test(int linear_solver, const Eigen::VectorXd y, const Eigen::ArrayXd weights, const Rcpp::NumericVector x, double rho, const Eigen::VectorXd adj_mean);
RcppExport SEXP _trendfilter_linear_single_solve_test(SEXP linear_solverSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP xSEXP, SEXP rhoSEXP, SEXP adj_meanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type linear_solver(linear_solverSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type adj_mean(adj_meanSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_single_solve_test(linear_solver, y, weights, x, rho, adj_mean));
    return rcpp_result_gen;
END_RCPP
}
// admm_lambda_seq
Rcpp::List admm_lambda_seq(NumericVector x, Eigen::VectorXd y, Eigen::ArrayXd weights, int k, Eigen::VectorXd lambda, int nlambda, double lambda_max, double lambda_min, double lambda_min_ratio, int max_iter, double rho_scale, double tol, int linear_solver, double space_tolerance_ratio);
RcppExport SEXP _trendfilter_admm_lambda_seq(SEXP xSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP kSEXP, SEXP lambdaSEXP, SEXP nlambdaSEXP, SEXP lambda_maxSEXP, SEXP lambda_minSEXP, SEXP lambda_min_ratioSEXP, SEXP max_iterSEXP, SEXP rho_scaleSEXP, SEXP tolSEXP, SEXP linear_solverSEXP, SEXP space_tolerance_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type y(ySEXP);
    Rcpp::traits::input_parameter< Eigen::ArrayXd >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type nlambda(nlambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min(lambda_minSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type rho_scale(rho_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type linear_solver(linear_solverSEXP);
    Rcpp::traits::input_parameter< double >::type space_tolerance_ratio(space_tolerance_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(admm_lambda_seq(x, y, weights, k, lambda, nlambda, lambda_max, lambda_min, lambda_min_ratio, max_iter, rho_scale, tol, linear_solver, space_tolerance_ratio));
    return rcpp_result_gen;
END_RCPP
}
// get_dk_mat
Eigen::SparseMatrix<double> get_dk_mat(int k, NumericVector xd, bool tf_weighting);
RcppExport SEXP _trendfilter_get_dk_mat(SEXP kSEXP, SEXP xdSEXP, SEXP tf_weightingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xd(xdSEXP);
    Rcpp::traits::input_parameter< bool >::type tf_weighting(tf_weightingSEXP);
    rcpp_result_gen = Rcpp::wrap(get_dk_mat(k, xd, tf_weighting));
    return rcpp_result_gen;
END_RCPP
}
// get_lambda_max
double get_lambda_max(const NumericVector& x, const Eigen::VectorXd& y, const Eigen::ArrayXd& weights, int k);
RcppExport SEXP _trendfilter_get_lambda_max(SEXP xSEXP, SEXP ySEXP, SEXP weightsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXd& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda_max(x, y, weights, k));
    return rcpp_result_gen;
END_RCPP
}
// calc_degrees_of_freedom
int calc_degrees_of_freedom(Eigen::VectorXd const& v, int k, double tol);
RcppExport SEXP _trendfilter_calc_degrees_of_freedom(SEXP vSEXP, SEXP kSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd const& >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_degrees_of_freedom(v, k, tol));
    return rcpp_result_gen;
END_RCPP
}
// get_lambda_seq_r
Eigen::VectorXd get_lambda_seq_r(Eigen::VectorXd lambda, double lambda_max, double lambda_min, double lambda_min_ratio, int n_lambda);
RcppExport SEXP _trendfilter_get_lambda_seq_r(SEXP lambdaSEXP, SEXP lambda_maxSEXP, SEXP lambda_minSEXP, SEXP lambda_min_ratioSEXP, SEXP n_lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_max(lambda_maxSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min(lambda_minSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min_ratio(lambda_min_ratioSEXP);
    Rcpp::traits::input_parameter< int >::type n_lambda(n_lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_lambda_seq_r(lambda, lambda_max, lambda_min, lambda_min_ratio, n_lambda));
    return rcpp_result_gen;
END_RCPP
}
// Dkv
Eigen::VectorXd Dkv(Eigen::VectorXd v, int k, const NumericVector& xd, bool tf_weighting);
RcppExport SEXP _trendfilter_Dkv(SEXP vSEXP, SEXP kSEXP, SEXP xdSEXP, SEXP tf_weightingSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type xd(xdSEXP);
    Rcpp::traits::input_parameter< bool >::type tf_weighting(tf_weightingSEXP);
    rcpp_result_gen = Rcpp::wrap(Dkv(v, k, xd, tf_weighting));
    return rcpp_result_gen;
END_RCPP
}
// is_equal_space
bool is_equal_space(Rcpp::NumericVector x, double space_tolerance_ratio);
RcppExport SEXP _trendfilter_is_equal_space(SEXP xSEXP, SEXP space_tolerance_ratioSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type space_tolerance_ratio(space_tolerance_ratioSEXP);
    rcpp_result_gen = Rcpp::wrap(is_equal_space(x, space_tolerance_ratio));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_trendfilter_smat_to_mat", (DL_FUNC) &_trendfilter_smat_to_mat, 3},
    {"_trendfilter_configure_denseD_test", (DL_FUNC) &_trendfilter_configure_denseD_test, 2},
    {"_trendfilter_linear_single_solve_test", (DL_FUNC) &_trendfilter_linear_single_solve_test, 6},
    {"_trendfilter_admm_lambda_seq", (DL_FUNC) &_trendfilter_admm_lambda_seq, 14},
    {"_trendfilter_get_dk_mat", (DL_FUNC) &_trendfilter_get_dk_mat, 3},
    {"_trendfilter_get_lambda_max", (DL_FUNC) &_trendfilter_get_lambda_max, 4},
    {"_trendfilter_calc_degrees_of_freedom", (DL_FUNC) &_trendfilter_calc_degrees_of_freedom, 3},
    {"_trendfilter_get_lambda_seq_r", (DL_FUNC) &_trendfilter_get_lambda_seq_r, 5},
    {"_trendfilter_Dkv", (DL_FUNC) &_trendfilter_Dkv, 4},
    {"_trendfilter_is_equal_space", (DL_FUNC) &_trendfilter_is_equal_space, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_trendfilter(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
