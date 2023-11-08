#ifndef __FFBS_H__
#define __FFBS_H__

#include "03_ffbs_helper.h"

Rcpp::List kf_ff(const arma::mat& yObs,
                 const arma::mat& uReg,
                 const arma::mat& wReg,
                 const int dimX,
                 const int dimY,
                 const int TT,
                 const arma::colvec& x00, 
                 const arma::mat& P00,
                 // u00,
                 const arma::mat& A, 
                 // const arma::mat& B,
                 const arma::mat& C,
                 const arma::mat& D, 
                 const arma::mat& Q,
                 const arma::mat& R,
                 bool PDSTORE = false);
arma::mat compute_mat_reg(const arma::mat& mat, const arma::mat& reg);
#endif