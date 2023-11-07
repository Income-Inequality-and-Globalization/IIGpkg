#include "01_gibbs_sampler.h"
//' Sum in the inverse of the posterior variance for loadings/partial effects
//'
//' @param availableObs available observations
//' @param id_f `integer` (sequence); index for selection of correct factors
//' @param fPost backward sampled states (FFBS output)
//' @param w_regs 
//' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
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
  // if (is.null(w_regs)) {
  //   for (tt in availableObs) {
  //     f <- fPost[, tt]
  //     V <- Viarray[, , tt]
  //     summ <- summ + kronecker(f %*% t(f), solve(V))
  //   }
  // } else {
  //   for (tt in availableObs) {
  //     f <- c(fPost[, tt], w_regs[, tt])
  //     V <- Viarray[, , tt]
  //     summ <- summ + kronecker(f %*% t(f), solve(V))
  //   }
  // // }
  // return(summ);
}
