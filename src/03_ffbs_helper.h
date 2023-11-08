#ifndef __FFBS_HELPER_H__
#define __FFBS_HELPER_H__

#include <RcppArmadillo.h>

arma::mat compute_mat_reg(const arma::mat& mat, const arma::mat& reg);
arma::mat copute_Xtt_1(const arma::mat& A, const arma::mat& xtt);
arma::mat compute_ptt_1(const arma::mat& A,
                        const arma::mat&Ptt,
                        const arma::mat&Q);
arma::mat compute_Lt(const arma::mat& C,
                     const arma::mat& Ptt1,
                     const arma::mat& R);
arma::mat compute_Kt(const arma::mat& Ptt1, const arma::mat& C);

#endif