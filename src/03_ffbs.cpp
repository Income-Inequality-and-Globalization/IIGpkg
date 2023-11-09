#include "03_ffbs.h"
//[[Rcpp::export]]
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
                 bool PDSTORE) {
  arma::mat xtt1STORE(dimX, TT + 1);
  arma::cube Ptt1STORE(dimX, dimX, TT + 1);
  bool AllObsMissing;
  arma::mat CAdj;
  arma::mat RAdj;
  arma::colvec yObsAdj;
  arma::colvec DwRegAdj;
  arma::mat C_adj;
  arma::mat R_adj;
  arma::mat yObs_adj;
  arma::mat DwReg_adj;
  arma::uvec MissingObs;
  arma::mat W;
  
  arma::mat xtt(dimX, TT);
  arma::cube Ptt(dimX, dimX, TT);
  arma::mat Lt, Kt;
  arma::vec kGain;
  if (PDSTORE) {
    // xtt1STORE <- matrix(0, nrow = dimX, ncol = TT + 1)
    // arma::mat xtt1STORE(dimX, TT + 1);
    // Ptt1STORE <- array(0, dim = c(dimX, dimX, TT + 1),
    //                    dimnames = list(NULL, NULL, as.character(1:(TT + 1))))
    // arma::cube Ptt1STORE(dimX, dimX, TT + 1);
  }
  arma::mat D_wreg = compute_mat_reg(D, wReg);
  // 
  // 
  //     xtt1 <- computeXtt1(A, x00, BuRegInit[, 1], dimX)
  arma::mat xtt1 = compute_Xtt_1(A, x00);
  //     Ptt1 <- computePtt1(A, P00, Q)
  arma::mat Ptt1 = compute_Ptt_1(A, P00, Q);
  // 
  if (PDSTORE) {
    xtt1STORE.col(0)   = xtt1;
    Ptt1STORE.slice(0) = Ptt1; 
  }
  
  for (arma::uword t = 0; t <  TT; ++t) {
  // arma::uword t = 0;
    AllObsMissing = false;
    CAdj = C;

    RAdj = R.slice(t);
    yObsAdj = yObs.col(t);
    DwRegAdj = D_wreg.col(t);
    MissingObs = arma::find_nan(yObsAdj);
    
    if (MissingObs.size() != 0) {
      if (MissingObs.size() == dimY) {
        AllObsMissing = true;
      } else {
        W.eye(dimY, dimY);
        W.shed_rows(MissingObs);
        
        for(arma::uword dy = 0; dy < dimY; ++dy) {
          for (auto mo : MissingObs) {
            CAdj.submat(mo, dy, mo, dy) = 0;
          }
        }
        for (auto mo2 : MissingObs) {
          yObsAdj(mo2) = 0;
          for (auto mo1 : MissingObs) {
            RAdj.submat(mo1, mo2, mo1, mo2) = 0;
          }
        }
        
        C_adj = W * CAdj;
        R_adj = W * RAdj * W.t();
        yObs_adj = W * yObsAdj;
        DwReg_adj = W * DwRegAdj;
      }
    } else {
      C_adj = CAdj;
      R_adj = RAdj;
      yObs_adj = yObsAdj;
      DwReg_adj = DwRegAdj;
    }
    if(!AllObsMissing) {
      Lt = compute_Lt(C_adj, Ptt1, R_adj);
      Kt = compute_Kt(Ptt1, C_adj);
      // period t quantities for current iteration
      kGain      = compute_kG(yObs_adj, C_adj, xtt1, DwReg_adj);
      xtt.col(t)   = compute_Xtt(xtt1, Kt, Lt, kGain);
      Ptt.slice(t) = compute_Ptt(Ptt1, Kt, Lt, C_adj, R_adj);
    } else {
      xtt.col(t) = xtt1;
      Ptt.slice(t) = Ptt1;
    }


    xtt1 = compute_Xtt_1(A, xtt.col(t));
    Ptt1 = compute_Ptt_1(A, Ptt.slice(t), Q);

    // period t+1 quantities for next iteration
    if (PDSTORE) {
      xtt1STORE.col(t + 1)   = xtt1;
      Ptt1STORE.slice(t + 1) = Ptt1;
    }
  }
  if (PDSTORE) {
    return(Rcpp::List::create(Rcpp::Named("mfdEXP") = xtt,
                              Rcpp::Named("mfdVAR") = Ptt,
                              Rcpp::Named("pddEXP") = xtt1STORE,
                              Rcpp::Named("pddVAR") = Ptt1STORE));
    
  } else {
    // return(Rcpp::List::create(Rcpp::Named("Lt") = Lt,
    //                           Rcpp::Named("mfdVAR") = Ptt));
    return(Rcpp::List::create(Rcpp::Named("mfdEXP") = xtt,
                              Rcpp::Named("mfdVAR") = Ptt));
  }
}
//' Backward Sampling based Kalman-Filter Output
//'
//' @param TT time
//' @param nfac Number of factors
//' @param Phi phi param
//' @param Q VCM part
//' @param filt_f Means of filtering distribution from Kalman-Filter
//' @param filt_P Variances of filtering distribution from Kalman-Filter
//'
//' @return posterior filtered states
//[[Rcpp::export]]
arma::mat bs(const int TT,
             const int nfac,
             const arma::mat& Phi,
             const arma::mat& Q,
             const arma::mat& filt_f,
             const arma::cube& filt_P) {
  arma::mat fPost(nfac, TT);

  arma::colvec fTT = filt_f.col(TT - 1);
  arma::colvec ftt(TT);

  arma::mat PTT = filt_P.slice(TT - 1);
  arma::mat Ptt(arma::size(PTT));

  fPost.col(TT - 1) = arma::mvnrnd(fTT, PTT);
  
  arma::mat IS(TT, TT);
  arma::mat sim_v(TT, TT);
  arma::colvec sim_m(TT);

  for (int t = TT - 2;  t >= 0; --t) {
    ftt = filt_f.col(t);
    Ptt = filt_P.slice(t);
    IS = Ptt * Phi.t() * arma::inv(Phi * Ptt * Phi.t() + Q);
    sim_m = ftt + IS * (fPost.col(t + 1) - Phi * ftt);
    sim_v = Ptt - IS * Phi * Ptt;
    sim_v = 0.5 * (sim_v + sim_v.t());
  // if(!matrixcalc::is.positive.definite(sim_v)){
  //   print(iter)
  //   print(sim_v)
  //   print(filt_P)
  // }
    // fPost[, t] = MASS::mvrnorm(n = 1, mu = sim_m, Sigma = sim_v)
    fPost.col(t) = arma::mvnrnd(sim_m, sim_v);
  }
  return(fPost);
//   
//   fPost = matrix(rep(0, nfac * TT), ncol = TT)
//   fPost[, TT] = MASS::mvrnorm(n = 1, mu = fTT, Sigma = PTT)
//   
//   for (t in (TT - 1):1) {
//     ftt = filt_f[, t]
//     Ptt = filt_P[, , t]
//     IS = Ptt %*% t(Phi) %*% solve(Phi %*% Ptt %*% t(Phi) + Q)
//     sim_m = ftt + IS %*% (fPost[, t + 1] - Phi %*% ftt)
//     sim_v = Ptt - IS %*% Phi %*% Ptt
//     sim_v = 0.5 * (sim_v + t(sim_v))
// // if(!matrixcalc::is.positive.definite(sim_v)){
// //   print(iter)
// //   print(sim_v)
// //   print(filt_P)
// // }
//     fPost[, t] = MASS::mvrnorm(n = 1, mu = sim_m, Sigma = sim_v)
//   }
//   
  // return(fPost);
}