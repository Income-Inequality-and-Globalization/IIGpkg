#include "01_gibbs_sampler.h"
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @param invOmega0 inverse prior matrix
//' @param availableObs available observations
//' @param id_f `integer` (sequence); index for selection of correct factors
//' @param fPost backward sampled states (FFBS output)
//' @param w_regs regressor matrix
//' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
//'
//' @return summation of kronecker products of appropriate dimension
//'
//[[Rcpp::export]]
arma::mat compute_B_mean_cpp(const arma::mat Omega,
                             const arma::colvec& invOmega_B0_D0,
                             const Rcpp::IntegerVector& availableObs,
                             const arma::mat selectR,
                             const arma::mat& fPost,
                             const arma::mat& w_regs,
                             const arma::mat yiObs,
                             const arma::cube& Viarray) {

  int f_size = fPost.n_rows + w_regs.n_rows;
  int v_size = Viarray.n_rows;
  // int w_size = arma::sum(arma::sum(w_regs));
  
  arma::colvec out_mean;

  arma::colvec out_rhs(f_size * v_size);
  out_rhs = invOmega_B0_D0 + sum_f_y_v(availableObs,
                                       fPost,
                                       w_regs,
                                       yiObs,
                                       Viarray);

  out_mean = Omega * (selectR * out_rhs);
  return(out_mean);
}
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @param invOmega0 inverse prior matrix
//' @param availableObs available observations
//' @param id_f `integer` (sequence); index for selection of correct factors
//' @param fPost backward sampled states (FFBS output)
//' @param w_regs regressor matrix
//' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
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
  out_omega += sum_ff_kron_v(availableObs, fPost, w_regs, Viarray);
  return(arma::inv(selectR * out_omega * selectR.t()));
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
//' @inheritParams compute_Omega1_cpp
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