#ifndef GIBBS_SAMPLER_H
#define GIBBS_SAMPLER_H

#include <RcppArmadillo.h>
#include "00_ghk_sampler.h"

arma::rowvec sample_B_D_cpp(const Rcpp::IntegerVector& availableObs,
                            const arma::mat& invOmega0,
                            const arma::colvec& B0_D0,
                            const arma::mat& fPost,
                            const arma::mat& w_regs,
                            const arma::cube& Viarray,
                            const arma::mat& yiObs,
                            const arma::mat& selectR,
                            const arma::rowvec& lower, 
                            const arma::rowvec& upper);
Rcpp::List compute_B_post_full_cpp(const Rcpp::IntegerVector& availableObs,
                                   const arma::mat& invOmega0,
                                   const arma::colvec& B0_D0,
                                   const arma::mat& fPost,
                                   const arma::mat& w_regs,
                                   const arma::cube& Viarray,
                                   const arma::mat& yiObs,
                                   const arma::mat& selectR);
arma::rowvec compute_B_mean_cpp(const arma::mat& Omega,
                                const arma::mat& invOmega0,
                                const arma::colvec& B0_D0,
                                const Rcpp::IntegerVector& availableObs,
                                const arma::mat& selectR,
                                const arma::mat& fPost,
                                const arma::mat& w_regs,
                                const arma::mat& yiObs,
                                const arma::cube& Viarray);
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
arma::mat compute_sigma_adjust_cpp(const arma::mat& Omega);
Rcpp::List get_BD_sample(const arma::rowvec& draw,
                         int num_jnt_fac, 
                         int num_y);

#endif