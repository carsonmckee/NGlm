// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_chains
Eigen::MatrixXd get_chains(const Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const Eigen::MatrixXd& X_t_X, const Eigen::MatrixXd& X_t_y, const int& n_samp, const int& n_thin, const int& n, const int& p, const double& M, const double& init_alpha, const Eigen::MatrixXd& init_beta, const Eigen::ArrayXd& init_psi, const double& init_sigma, const double& init_lambda, const double& init_gamma, const double& init_prop_sd, const int& tuning_time, const int& tuning_freq, const int& verbose);
RcppExport SEXP _NGlm_get_chains(SEXP ySEXP, SEXP xSEXP, SEXP X_t_XSEXP, SEXP X_t_ySEXP, SEXP n_sampSEXP, SEXP n_thinSEXP, SEXP nSEXP, SEXP pSEXP, SEXP MSEXP, SEXP init_alphaSEXP, SEXP init_betaSEXP, SEXP init_psiSEXP, SEXP init_sigmaSEXP, SEXP init_lambdaSEXP, SEXP init_gammaSEXP, SEXP init_prop_sdSEXP, SEXP tuning_timeSEXP, SEXP tuning_freqSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X_t_X(X_t_XSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type X_t_y(X_t_ySEXP);
    Rcpp::traits::input_parameter< const int& >::type n_samp(n_sampSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_thin(n_thinSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const double& >::type M(MSEXP);
    Rcpp::traits::input_parameter< const double& >::type init_alpha(init_alphaSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type init_beta(init_betaSEXP);
    Rcpp::traits::input_parameter< const Eigen::ArrayXd& >::type init_psi(init_psiSEXP);
    Rcpp::traits::input_parameter< const double& >::type init_sigma(init_sigmaSEXP);
    Rcpp::traits::input_parameter< const double& >::type init_lambda(init_lambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type init_gamma(init_gammaSEXP);
    Rcpp::traits::input_parameter< const double& >::type init_prop_sd(init_prop_sdSEXP);
    Rcpp::traits::input_parameter< const int& >::type tuning_time(tuning_timeSEXP);
    Rcpp::traits::input_parameter< const int& >::type tuning_freq(tuning_freqSEXP);
    Rcpp::traits::input_parameter< const int& >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(get_chains(y, x, X_t_X, X_t_y, n_samp, n_thin, n, p, M, init_alpha, init_beta, init_psi, init_sigma, init_lambda, init_gamma, init_prop_sd, tuning_time, tuning_freq, verbose));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_NGlm_get_chains", (DL_FUNC) &_NGlm_get_chains, 19},
    {NULL, NULL, 0}
};

RcppExport void R_init_NGlm(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}