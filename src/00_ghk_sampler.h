#ifndef __GHK_SAMPLER_H__
#define __GHK_SAMPLER_H__

#include "00_ghk_sampler_helper.h"

arma::mat rtmvnormcpp(const arma::mat& mean, 
                      const arma::mat& sigma, 
                      const arma::mat& blc,
                      const arma::mat& lower, 
                      const arma::mat& upper,
                      const arma::mat& init,
                      const arma::uword burn = 10);

#endif