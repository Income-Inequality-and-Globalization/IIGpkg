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
// sample_B_D_cpp
Rcpp::List sample_B_D_cpp(const Rcpp::IntegerVector& availableObs, const arma::mat& invOmega0, const arma::colvec& B0_D0, const arma::mat& fPost, const arma::mat& w_regs, const arma::cube& Viarray, const arma::mat& yiObs, const arma::mat& selectR, const arma::rowvec& lower, const arma::rowvec& upper, const int num_jnt_fac, const int num_y);
RcppExport SEXP _IIGpkg_sample_B_D_cpp(SEXP availableObsSEXP, SEXP invOmega0SEXP, SEXP B0_D0SEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP ViarraySEXP, SEXP yiObsSEXP, SEXP selectRSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP num_jnt_facSEXP, SEXP num_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invOmega0(invOmega0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type B0_D0(B0_D0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Viarray(ViarraySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yiObs(yiObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type selectR(selectRSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const int >::type num_jnt_fac(num_jnt_facSEXP);
    Rcpp::traits::input_parameter< const int >::type num_y(num_ySEXP);
    rcpp_result_gen = Rcpp::wrap(sample_B_D_cpp(availableObs, invOmega0, B0_D0, fPost, w_regs, Viarray, yiObs, selectR, lower, upper, num_jnt_fac, num_y));
    return rcpp_result_gen;
END_RCPP
}
// get_BD_sample
Rcpp::List get_BD_sample(const arma::rowvec& draw, int num_jnt_fac, int num_y);
RcppExport SEXP _IIGpkg_get_BD_sample(SEXP drawSEXP, SEXP num_jnt_facSEXP, SEXP num_ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::rowvec& >::type draw(drawSEXP);
    Rcpp::traits::input_parameter< int >::type num_jnt_fac(num_jnt_facSEXP);
    Rcpp::traits::input_parameter< int >::type num_y(num_ySEXP);
    rcpp_result_gen = Rcpp::wrap(get_BD_sample(draw, num_jnt_fac, num_y));
    return rcpp_result_gen;
END_RCPP
}
// compute_B_post_full_cpp
Rcpp::List compute_B_post_full_cpp(const Rcpp::IntegerVector& availableObs, const arma::mat& invOmega0, const arma::colvec& B0_D0, const arma::mat& fPost, const arma::mat& w_regs, const arma::cube& Viarray, const arma::mat& yiObs, const arma::mat& selectR);
RcppExport SEXP _IIGpkg_compute_B_post_full_cpp(SEXP availableObsSEXP, SEXP invOmega0SEXP, SEXP B0_D0SEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP ViarraySEXP, SEXP yiObsSEXP, SEXP selectRSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invOmega0(invOmega0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type B0_D0(B0_D0SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Viarray(ViarraySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yiObs(yiObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type selectR(selectRSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_B_post_full_cpp(availableObs, invOmega0, B0_D0, fPost, w_regs, Viarray, yiObs, selectR));
    return rcpp_result_gen;
END_RCPP
}
// compute_B_mean_cpp
arma::rowvec compute_B_mean_cpp(const arma::mat& Omega, const arma::mat& invOmega0, const arma::colvec& B0_D0, const Rcpp::IntegerVector& availableObs, const arma::mat& selectR, const arma::mat& fPost, const arma::mat& w_regs, const arma::mat& yiObs, const arma::cube& Viarray);
RcppExport SEXP _IIGpkg_compute_B_mean_cpp(SEXP OmegaSEXP, SEXP invOmega0SEXP, SEXP B0_D0SEXP, SEXP availableObsSEXP, SEXP selectRSEXP, SEXP fPostSEXP, SEXP w_regsSEXP, SEXP yiObsSEXP, SEXP ViarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type invOmega0(invOmega0SEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type B0_D0(B0_D0SEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type availableObs(availableObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type selectR(selectRSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fPost(fPostSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w_regs(w_regsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type yiObs(yiObsSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type Viarray(ViarraySEXP);
    rcpp_result_gen = Rcpp::wrap(compute_B_mean_cpp(Omega, invOmega0, B0_D0, availableObs, selectR, fPost, w_regs, yiObs, Viarray));
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
// ffbs
arma::mat ffbs(const arma::mat& yObs, const arma::mat& wReg, const arma::uword dimX, const arma::uword dimY, const arma::uword TT, const arma::colvec& x00, const arma::mat& P00, const arma::mat& A, const arma::mat& C, const arma::mat& D, const arma::mat& Q, const arma::cube& R_reordered, bool PDSTORE);
RcppExport SEXP _IIGpkg_ffbs(SEXP yObsSEXP, SEXP wRegSEXP, SEXP dimXSEXP, SEXP dimYSEXP, SEXP TTSEXP, SEXP x00SEXP, SEXP P00SEXP, SEXP ASEXP, SEXP CSEXP, SEXP DSEXP, SEXP QSEXP, SEXP R_reorderedSEXP, SEXP PDSTORESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type yObs(yObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wReg(wRegSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type dimX(dimXSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type dimY(dimYSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type x00(x00SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P00(P00SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type R_reordered(R_reorderedSEXP);
    Rcpp::traits::input_parameter< bool >::type PDSTORE(PDSTORESEXP);
    rcpp_result_gen = Rcpp::wrap(ffbs(yObs, wReg, dimX, dimY, TT, x00, P00, A, C, D, Q, R_reordered, PDSTORE));
    return rcpp_result_gen;
END_RCPP
}
// kf_ff
Rcpp::List kf_ff(const arma::mat& yObs, const arma::mat& wReg, const arma::uword dimX, const arma::uword dimY, const arma::uword TT, const arma::colvec& x00, const arma::mat& P00, const arma::mat& A, const arma::mat& C, const arma::mat& D, const arma::mat& Q, const arma::cube& R, bool PDSTORE, bool LLVALUE);
RcppExport SEXP _IIGpkg_kf_ff(SEXP yObsSEXP, SEXP wRegSEXP, SEXP dimXSEXP, SEXP dimYSEXP, SEXP TTSEXP, SEXP x00SEXP, SEXP P00SEXP, SEXP ASEXP, SEXP CSEXP, SEXP DSEXP, SEXP QSEXP, SEXP RSEXP, SEXP PDSTORESEXP, SEXP LLVALUESEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type yObs(yObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type wReg(wRegSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type dimX(dimXSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type dimY(dimYSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const arma::colvec& >::type x00(x00SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type P00(P00SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type R(RSEXP);
    Rcpp::traits::input_parameter< bool >::type PDSTORE(PDSTORESEXP);
    Rcpp::traits::input_parameter< bool >::type LLVALUE(LLVALUESEXP);
    rcpp_result_gen = Rcpp::wrap(kf_ff(yObs, wReg, dimX, dimY, TT, x00, P00, A, C, D, Q, R, PDSTORE, LLVALUE));
    return rcpp_result_gen;
END_RCPP
}
// bs
arma::mat bs(const arma::uword TT, const int nfac, const arma::mat& Phi, const arma::mat& Q, const arma::mat& filt_f, const arma::cube& filt_P);
RcppExport SEXP _IIGpkg_bs(SEXP TTSEXP, SEXP nfacSEXP, SEXP PhiSEXP, SEXP QSEXP, SEXP filt_fSEXP, SEXP filt_PSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::uword >::type TT(TTSEXP);
    Rcpp::traits::input_parameter< const int >::type nfac(nfacSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type filt_f(filt_fSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type filt_P(filt_PSEXP);
    rcpp_result_gen = Rcpp::wrap(bs(TT, nfac, Phi, Q, filt_f, filt_P));
    return rcpp_result_gen;
END_RCPP
}
// compute_mat_reg
arma::mat compute_mat_reg(const arma::mat& mat, const arma::mat& reg);
RcppExport SEXP _IIGpkg_compute_mat_reg(SEXP matSEXP, SEXP regSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type reg(regSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_mat_reg(mat, reg));
    return rcpp_result_gen;
END_RCPP
}
// compute_Xtt_1
arma::mat compute_Xtt_1(const arma::mat& A, const arma::mat& xtt);
RcppExport SEXP _IIGpkg_compute_Xtt_1(SEXP ASEXP, SEXP xttSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xtt(xttSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Xtt_1(A, xtt));
    return rcpp_result_gen;
END_RCPP
}
// compute_Ptt_1
arma::mat compute_Ptt_1(const arma::mat& A, const arma::mat& Ptt, const arma::mat& Q);
RcppExport SEXP _IIGpkg_compute_Ptt_1(SEXP ASEXP, SEXP PttSEXP, SEXP QSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Ptt(PttSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Q(QSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Ptt_1(A, Ptt, Q));
    return rcpp_result_gen;
END_RCPP
}
// compute_Lt
arma::mat compute_Lt(const arma::mat& C, const arma::mat& Ptt1, const arma::mat& R);
RcppExport SEXP _IIGpkg_compute_Lt(SEXP CSEXP, SEXP Ptt1SEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Ptt1(Ptt1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Lt(C, Ptt1, R));
    return rcpp_result_gen;
END_RCPP
}
// compute_Kt
arma::mat compute_Kt(const arma::mat& Ptt1, const arma::mat& C);
RcppExport SEXP _IIGpkg_compute_Kt(SEXP Ptt1SEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Ptt1(Ptt1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Kt(Ptt1, C));
    return rcpp_result_gen;
END_RCPP
}
// compute_kG
arma::mat compute_kG(const arma::mat& yObs, const arma::mat& C, const arma::mat& xtt1, const arma::mat& DwReg);
RcppExport SEXP _IIGpkg_compute_kG(SEXP yObsSEXP, SEXP CSEXP, SEXP xtt1SEXP, SEXP DwRegSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type yObs(yObsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type xtt1(xtt1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DwReg(DwRegSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_kG(yObs, C, xtt1, DwReg));
    return rcpp_result_gen;
END_RCPP
}
// compute_Xtt
arma::mat compute_Xtt(const arma::mat& xtt1, const arma::mat& Kt, const arma::mat& Lt, const arma::mat& kGain);
RcppExport SEXP _IIGpkg_compute_Xtt(SEXP xtt1SEXP, SEXP KtSEXP, SEXP LtSEXP, SEXP kGainSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type xtt1(xtt1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Kt(KtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lt(LtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type kGain(kGainSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Xtt(xtt1, Kt, Lt, kGain));
    return rcpp_result_gen;
END_RCPP
}
// compute_Ptt
arma::mat compute_Ptt(const arma::mat& Ptt1, const arma::mat& Kt, const arma::mat& Lt, const arma::mat& C, const arma::mat& R);
RcppExport SEXP _IIGpkg_compute_Ptt(SEXP Ptt1SEXP, SEXP KtSEXP, SEXP LtSEXP, SEXP CSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Ptt1(Ptt1SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Kt(KtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lt(LtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type C(CSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Ptt(Ptt1, Kt, Lt, C, R));
    return rcpp_result_gen;
END_RCPP
}
// compute_b_diag_by_time
arma::cube compute_b_diag_by_time(const arma::cube& V_hat_array_A, const arma::uword num_y, const arma::uword NN, const arma::uword TT);
RcppExport SEXP _IIGpkg_compute_b_diag_by_time(SEXP V_hat_array_ASEXP, SEXP num_ySEXP, SEXP NNSEXP, SEXP TTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::cube& >::type V_hat_array_A(V_hat_array_ASEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type num_y(num_ySEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type NN(NNSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type TT(TTSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_b_diag_by_time(V_hat_array_A, num_y, NN, TT));
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
    {"_IIGpkg_sample_B_D_cpp", (DL_FUNC) &_IIGpkg_sample_B_D_cpp, 12},
    {"_IIGpkg_get_BD_sample", (DL_FUNC) &_IIGpkg_get_BD_sample, 3},
    {"_IIGpkg_compute_B_post_full_cpp", (DL_FUNC) &_IIGpkg_compute_B_post_full_cpp, 8},
    {"_IIGpkg_compute_B_mean_cpp", (DL_FUNC) &_IIGpkg_compute_B_mean_cpp, 9},
    {"_IIGpkg_compute_Omega1_cpp", (DL_FUNC) &_IIGpkg_compute_Omega1_cpp, 6},
    {"_IIGpkg_sum_ff_kron_v", (DL_FUNC) &_IIGpkg_sum_ff_kron_v, 4},
    {"_IIGpkg_sum_f_y_v", (DL_FUNC) &_IIGpkg_sum_f_y_v, 5},
    {"_IIGpkg_ffbs", (DL_FUNC) &_IIGpkg_ffbs, 13},
    {"_IIGpkg_kf_ff", (DL_FUNC) &_IIGpkg_kf_ff, 14},
    {"_IIGpkg_bs", (DL_FUNC) &_IIGpkg_bs, 6},
    {"_IIGpkg_compute_mat_reg", (DL_FUNC) &_IIGpkg_compute_mat_reg, 2},
    {"_IIGpkg_compute_Xtt_1", (DL_FUNC) &_IIGpkg_compute_Xtt_1, 2},
    {"_IIGpkg_compute_Ptt_1", (DL_FUNC) &_IIGpkg_compute_Ptt_1, 3},
    {"_IIGpkg_compute_Lt", (DL_FUNC) &_IIGpkg_compute_Lt, 3},
    {"_IIGpkg_compute_Kt", (DL_FUNC) &_IIGpkg_compute_Kt, 2},
    {"_IIGpkg_compute_kG", (DL_FUNC) &_IIGpkg_compute_kG, 4},
    {"_IIGpkg_compute_Xtt", (DL_FUNC) &_IIGpkg_compute_Xtt, 4},
    {"_IIGpkg_compute_Ptt", (DL_FUNC) &_IIGpkg_compute_Ptt, 5},
    {"_IIGpkg_compute_b_diag_by_time", (DL_FUNC) &_IIGpkg_compute_b_diag_by_time, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_IIGpkg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
