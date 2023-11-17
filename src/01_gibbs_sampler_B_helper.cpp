#include "01_gibbs_sampler.h"
//' Sample gkh-style (from truncated normal)
//'
//' @param availableObs available observations
//' @param invOmega0 inverse prior matrix
//' @param invOmega0_B0_D0 inverse prior matrix multiplied by B0,D0
//' @param fPost matrix of factors (backward sampled states via FFBS)
//' @param w_regs regressor matrix
//' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
//' @param yiObs matrix of obersevations
//' @param selectR selector matrix for missing y
//'
//' @return a list of two elements, the first being the posterior mean and the 
//'    second the VCM matrix for sampling
//'
//[[Rcpp::export]]
Rcpp::List sample_B_D_cpp(const Rcpp::IntegerVector& availableObs,
                          const arma::mat& invOmega0,
                          const arma::colvec& B0_D0,
                          const arma::mat& fPost,
                          const arma::mat& w_regs,
                          const arma::cube& Viarray,
                          const arma::mat& yiObs,
                          const arma::mat& selectR,
                          const arma::rowvec& lower, 
                          const arma::rowvec& upper,
                          const int num_jnt_fac,
                          const int num_y) {
  Rcpp::List out_sample_list;
  Rcpp::List out_mean_sigma;
  out_mean_sigma = compute_B_post_full_cpp(availableObs,
                                           invOmega0,
                                           B0_D0,
                                           fPost,
                                           w_regs,
                                           Viarray,
                                           yiObs,
                                           selectR);
  arma::rowvec mean_ghk = out_mean_sigma["B_mean"];
  arma::mat Sigma_ghk = out_mean_sigma["B_Sigma"];

  arma::mat tmp_identity = arma::eye(arma::size(Sigma_ghk));
  arma::rowvec tmp_zeros(mean_ghk.size());
  arma::rowvec tmp_sample(mean_ghk.size());
  
  // IIGpkg:::rtmvnormcpp(matrix(bmean, nrow = 1), Sigma, diag(ncol(Sigma)),
  //                      matrix(lower, nrow = 1), matrix(upper, nrow = 1),
  //                      init = matrix(0, ncol = ncol(Sigma)), burn = 10)
  tmp_sample = rtmvnormcpp(mean_ghk, Sigma_ghk, tmp_identity, lower, upper, tmp_zeros, 10);
  out_sample_list = get_BD_sample(tmp_sample, num_jnt_fac, num_y);
  return(out_sample_list);
}
//' Helper to separate joint, idiosyncratic and regressor samples
//'
//' @param draw full sample draw from truncated normal
//' @param num_jnt_fac an integer giving the number of joint factors
//' @param num_y integer giving the number of components in the measurements
//'
//' @return a list of thre elements, the first being the idiosyncratic factors,
//'    the second being the common factor part and the last being the regressors
//'
//[[Rcpp::export]]
Rcpp::List get_BD_sample(const arma::rowvec& draw,
                         int num_jnt_fac, 
                         int num_y) {
  Rcpp::List out;
  int separate_id_B_D = (num_jnt_fac + 1) * num_y;
  arma::rowvec Bvec = draw.subvec(0, separate_id_B_D - 1);
  arma::rowvec Dvec = draw.subvec(separate_id_B_D, draw.size() - 1);
  if (num_jnt_fac == 0) {
    out = Rcpp::List::create(
      Rcpp::Named("Bsamp_idi") = Bvec,
      Rcpp::Named("Bsamp_jnt") = R_NilValue,
      Rcpp::Named("Dsamp") = Dvec);
  } else if (num_jnt_fac == 1) {
    out = Rcpp::List::create(
      Rcpp::Named("Bsamp_idi") = Bvec.cols(num_jnt_fac * num_y, Bvec.size() - 1),
      Rcpp::Named("Bsamp_jnt") = Bvec.cols(0, num_jnt_fac * num_y - 1),
      Rcpp::Named("Dsamp") = Dvec);
  }
  return(out);
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @inheritParams sample_B_D_cpp
//'
//' @return a list of two elements, the first being the posterior mean and the 
//'    second the VCM matrix for sampling
//'
//[[Rcpp::export]]
Rcpp::List compute_B_post_full_cpp(const Rcpp::IntegerVector& availableObs,
                                   const arma::mat& invOmega0,
                                   const arma::colvec& B0_D0,
                                   const arma::mat& fPost,
                                   const arma::mat& w_regs,
                                   const arma::cube& Viarray,
                                   const arma::mat& yiObs,
                                   const arma::mat& selectR) {
  arma::mat B_Omega;
  arma::rowvec B_mean;
  B_Omega = compute_Omega1_cpp(invOmega0 ,
                               availableObs,
                               selectR,
                               fPost,
                               w_regs,
                               Viarray);
  B_mean = compute_B_mean_cpp(B_Omega,
                              invOmega0,
                              B0_D0,
                              availableObs, selectR,
                              fPost, w_regs,
                              yiObs,  Viarray);
  B_Omega = compute_sigma_adjust_cpp(B_Omega);
  
  return(Rcpp::List::create(Rcpp::Named("B_mean") = B_mean,
                            Rcpp::Named("B_Sigma") = B_Omega));
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @inheritParams sample_B_D_cpp
//'
//' @return summation of kronecker products of appropriate dimension
//'
//[[Rcpp::export]]
arma::rowvec compute_B_mean_cpp(const arma::mat& Omega,
                                const arma::mat& invOmega0,
                                const arma::colvec& B0_D0,
                                const Rcpp::IntegerVector& availableObs,
                                const arma::mat& selectR,
                                const arma::mat& fPost,
                                const arma::mat& w_regs,
                                const arma::mat& yiObs,
                                const arma::cube& Viarray) {
  int f_size = fPost.n_rows + w_regs.n_rows;
  int v_size = Viarray.n_rows;
  // int w_size = arma::sum(arma::sum(w_regs));
  arma::colvec out(f_size * v_size);
  out = invOmega0 * selectR * B0_D0;
  out += selectR * sum_f_y_v(availableObs, fPost, w_regs, yiObs, Viarray);
  out = Omega * out;
  return(out.t());
}
arma::mat compute_sigma_adjust_cpp(const arma::mat& Omega) {
  return(0.5 * Omega + 0.5 * Omega.t());
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @inheritParams sample_B_D_cpp
//'
//' @return summation of kronecker products of appropriate dimension
//'
//[[Rcpp::export]]
arma::mat compute_Omega1_cpp(const arma::mat& invOmega0,
                             const Rcpp::IntegerVector& availableObs,
                             const arma::mat selectR,
                             const arma::mat& fPost,
                             const arma::mat& w_regs,
                             const arma::cube& Viarray) {
  arma::mat out_omega = invOmega0;
  out_omega += selectR * sum_ff_kron_v(availableObs, fPost, w_regs, Viarray) * selectR.t();
  return(arma::inv(out_omega));
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @inheritParams compute_Omega1_cpp
//'
//' @return summation of kronecker products of appropriate dimension
//'
//[[Rcpp::export]]
arma::mat sum_ff_kron_v(const Rcpp::IntegerVector& availableObs,
                        const arma::mat fPost,
                        const arma::mat w_regs,
                        const arma::cube Viarray) {
  const Rcpp::IntegerVector tt_rng = availableObs - 1;
  int f_size = fPost.n_rows + w_regs.n_rows;
  int v_size = Viarray.n_rows;
  arma::mat f(f_size, tt_rng.length());
  f = arma::join_vert(fPost, w_regs);
  
  arma::mat out_ffkonv(f_size * v_size, f_size * v_size);
  for (auto tt : tt_rng) {
    out_ffkonv += arma::kron(f.col(tt) * f.col(tt).t(), Viarray.slice(tt).i());
  }
  return(out_ffkonv);
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @inheritParams sample_B_D_cpp
//'
//' @return summation of kronecker products of appropriate dimension
//'
//[[Rcpp::export]]
arma::colvec sum_f_y_v(const Rcpp::IntegerVector& availableObs,
                       const arma::mat fPost,
                       const arma::mat w_regs,
                       const arma::mat yiObs, 
                       const arma::cube Viarray) {


  const Rcpp::IntegerVector tt_rng = availableObs - 1;
  int f_size = fPost.n_rows + w_regs.n_rows;
  int v_size = Viarray.n_rows;
  arma::mat f(f_size, tt_rng.length());
  f = arma::join_vert(fPost, w_regs);
  arma::mat out_fyv(v_size, f_size);
  for (auto tt : tt_rng) {
    out_fyv += Viarray.slice(tt).i() * yiObs.col(tt) * f.col(tt).t();
  }
  return(out_fyv.reshape(v_size * f_size, 1));
}