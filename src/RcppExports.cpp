// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rtmvnormcpp
arma::mat rtmvnormcpp(const arma::mat& mean, const arma::mat& sigma, const arma::mat& blc, const arma::mat& lower, const arma::mat& upper, const arma::mat& init, const arma::uword burn);
RcppExport SEXP _IIGpkg_rtmvnormcpp(SEXP meanSEXP, SEXP sigmaSEXP, SEXP blcSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP initSEXP, SEXP burnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type blc(blcSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type init(initSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type burn(burnSEXP);
    rcpp_result_gen = Rcpp::wrap(rtmvnormcpp(mean, sigma, blc, lower, upper, init, burn));
    return rcpp_result_gen;
END_RCPP
}
// norm_rej
double norm_rej(const double a, const double b);
RcppExport SEXP _IIGpkg_norm_rej(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(norm_rej(a, b));
    return rcpp_result_gen;
END_RCPP
}
// unif_rej
double unif_rej(const double a, const double b);
RcppExport SEXP _IIGpkg_unif_rej(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(unif_rej(a, b));
    return rcpp_result_gen;
END_RCPP
}
// halfnorm_rej
double halfnorm_rej(const double a, const double b);
RcppExport SEXP _IIGpkg_halfnorm_rej(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(halfnorm_rej(a, b));
    return rcpp_result_gen;
END_RCPP
}
// exp_rej
double exp_rej(const double a, const double b);
RcppExport SEXP _IIGpkg_exp_rej(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_rej(a, b));
    return rcpp_result_gen;
END_RCPP
}
// rtnormcpp
arma::vec rtnormcpp(const arma::vec& mean, const double sd, const arma::vec& lower, const arma::vec& upper);
RcppExport SEXP _IIGpkg_rtnormcpp(SEXP meanSEXP, SEXP sdSEXP, SEXP lowerSEXP, SEXP upperSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< const double >::type sd(sdSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type upper(upperSEXP);
    rcpp_result_gen = Rcpp::wrap(rtnormcpp(mean, sd, lower, upper));
    return rcpp_result_gen;
END_RCPP
}
// compute_B_mean_cpp
arma::mat compute_B_mean_cpp(const arma::mat Omega, const arma::colvec& invOmega_B0_D0, const Rcpp::IntegerVector& availableObs, const arma::mat selectR, const arma::mat& fPost, const arma::mat& w_regs, const arma::mat yiObs, const arma::cube& Viarray);
RcppExport SEXP _IIGpkg_compute_B_mean_cpp(SEXP OmegaSEXP, SEXP invOmega_B0_D0SEXP, SEXP availableObsSEXP, SEXP selectRSEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP yiObsSEXP, SEXP ViarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type invOmega_B0_D0(invOmega_B0_D0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type selectR(selectRSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type yiObs(yiObsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Viarray(ViarraySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_B_mean_cpp(Omega, invOmega_B0_D0, availableObs, selectR, fPost, w_regs, yiObs, Viarray));
    return rcpp_result_gen;
END_RCPP
}
// compute_Omega1_cpp
arma::mat compute_Omega1_cpp(const arma::mat& invOmega0, const Rcpp::IntegerVector& availableObs, const arma::mat selectR, const arma::mat& fPost, const arma::mat& w_regs, const arma::cube& Viarray);
RcppExport SEXP _IIGpkg_compute_Omega1_cpp(SEXP invOmega0SEXP, SEXP availableObsSEXP, SEXP selectRSEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP ViarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type invOmega0(invOmega0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type selectR(selectRSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Viarray(ViarraySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Omega1_cpp(invOmega0, availableObs, selectR, fPost, w_regs, Viarray));
    return rcpp_result_gen;
END_RCPP
}
// sum_ff_kron_v
arma::mat sum_ff_kron_v(const Rcpp::IntegerVector& availableObs, const arma::mat fPost, const arma::mat w_regs, const arma::cube Viarray);
RcppExport SEXP _IIGpkg_sum_ff_kron_v(SEXP availableObsSEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP ViarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Viarray(ViarraySEXP);
    rcpp_result_gen = Rcpp::wrap(sum_ff_kron_v(availableObs, fPost, w_regs, Viarray));
    return rcpp_result_gen;
END_RCPP
}
// sum_f_y_v
arma::colvec sum_f_y_v(const Rcpp::IntegerVector& availableObs, const arma::mat fPost, const arma::mat w_regs, const arma::mat yiObs, const arma::cube Viarray);
RcppExport SEXP _IIGpkg_sum_f_y_v(SEXP availableObsSEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP yiObsSEXP, SEXP ViarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type yiObs(yiObsSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Viarray(ViarraySEXP);
    rcpp_result_gen = Rcpp::wrap(sum_f_y_v(availableObs, fPost, w_regs, yiObs, Viarray));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IIGpkg_rtmvnormcpp", (DL_FUNC) &_IIGpkg_rtmvnormcpp, 7},
    {"_IIGpkg_norm_rej", (DL_FUNC) &_IIGpkg_norm_rej, 2},
    {"_IIGpkg_unif_rej", (DL_FUNC) &_IIGpkg_unif_rej, 2},
    {"_IIGpkg_halfnorm_rej", (DL_FUNC) &_IIGpkg_halfnorm_rej, 2},
    {"_IIGpkg_exp_rej", (DL_FUNC) &_IIGpkg_exp_rej, 2},
    {"_IIGpkg_rtnormcpp", (DL_FUNC) &_IIGpkg_rtnormcpp, 4},
    {"_IIGpkg_compute_B_mean_cpp", (DL_FUNC) &_IIGpkg_compute_B_mean_cpp, 8},
    {"_IIGpkg_compute_Omega1_cpp", (DL_FUNC) &_IIGpkg_compute_Omega1_cpp, 6},
    {"_IIGpkg_sum_ff_kron_v", (DL_FUNC) &_IIGpkg_sum_ff_kron_v, 4},
    {"_IIGpkg_sum_f_y_v", (DL_FUNC) &_IIGpkg_sum_f_y_v, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_IIGpkg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
