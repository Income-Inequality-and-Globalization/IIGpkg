#ifndef GIBBS_SAMPLER_H
#define GIBBS_SAMPLER_H

#include <RcppArmadillo.h>

arma::mat sum_ff_kron_v(const Rcpp::IntegerVector& availableObs,
                        const arma::mat fPost,
                        const arma::mat w_regs,
                        const arma::cube Viarray);

#endif