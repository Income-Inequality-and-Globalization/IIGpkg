#ifndef GIBBS_SAMPLER_H
#define GIBBS_SAMPLER_H

#include <RcppArmadillo.h>

arma::mat compute_Omega1_cpp(const arma::mat& invOmega0,
                             const Rcpp::IntegerVector& availableObs,
                             const arma::mat selectR,
                             const arma::mat& fPost,
                             const arma::mat& w_regs,
                             const arma::cube& Viarray);
arma::mat sum_ff_kron_v(const Rcpp::IntegerVector& availableObs,
                        const arma::mat fPost,
                        const arma::mat w_regs,
                        const arma::cube Viarray);
arma::colvec sum_f_y_v(const Rcpp::IntegerVector& availableObs,
                       const arma::mat fPost,
                       const arma::mat w_regs,
                       const arma::mat yiObs, 
                       const arma::cube Viarray);

#endif