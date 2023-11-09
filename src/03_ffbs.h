#ifndef __FFBS_H__
#define __FFBS_H__

#include "03_ffbs_helper.h"

Rcpp::List kf_ff(const arma::mat& yObs,
                 const arma::mat& wReg,
                 const arma::uword dimX,
                 const arma::uword dimY,
                 const arma::uword TT,
                 const arma::colvec& x00, 
                 const arma::mat& P00,
                 const arma::mat& A, 
                 const arma::mat& C,
                 const arma::mat& D, 
                 const arma::mat& Q,
                 const arma::cube& R,
                 bool PDSTORE);
arma::mat bs(const int TT,
             const int nfac,
             const arma::mat& Phi,
             const arma::mat& Q,
             const arma::mat& filt_f,
             const arma::cube& filt_P);
#endif