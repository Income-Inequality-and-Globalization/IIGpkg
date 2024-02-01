#include "03_ffbs.h"

//[[Rcpp::export]]
arma::mat ffbs(const arma::mat& yObs,
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
               const arma::cube& R_reordered,
               bool PDSTORE) {
  const arma::uword num_y = R_reordered.n_rows;
  const arma::uword NN = dimY / R_reordered.n_rows;
  arma::cube R = compute_b_diag_by_time(R_reordered, num_y, NN, TT);

  arma::mat mfdEXP;
  arma::mat mfdVAR;
  Rcpp::List out_kf = Rcpp::List::create(Rcpp::Named("mfdEXP") = mfdEXP,
                                         Rcpp::Named("mfdVAR") = mfdVAR);
  out_kf = kf_ff(yObs, wReg, dimX, dimY, TT, x00, P00, A, C, D, Q, R, PDSTORE, false);
  return(bs(TT, dimX, A, Q, out_kf["mfdEXP"], out_kf["mfdVAR"]));
}
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
                 bool PDSTORE,
                 bool LLVALUE) {
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

  arma::mat D_wreg = compute_mat_reg(D, wReg);
  arma::mat xtt1 = compute_Xtt_1(A, x00);
  arma::mat Ptt1 = compute_Ptt_1(A, P00, Q);
  
  // double logDetVarY;
  double ll_out_tmp;
  double ll_out;
  // double dim_tmp;
  // double acc;
  double missing_obs_total;
  arma::uvec indices_nan_total = arma::find_nan(yObs);
  arma::vec vec_tmp_ll;
  arma::colvec meanY;
  arma::mat VarY;
  arma::mat chol_VarY;
  arma::vec chol_VarY_diag;
  // 
  if (PDSTORE) {
    xtt1STORE.col(0)   = xtt1;
    Ptt1STORE.slice(0) = Ptt1; 
  }
  //
  if (LLVALUE) {
    // logDetVarY = 0;
    ll_out_tmp = 0;
    ll_out = 0;
    // dim_tmp = 0;
    // acc = 0;
    missing_obs_total = indices_nan_total.size();
  }
  //
  for (arma::uword t = 0; t <  TT; ++t) {
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

    if(LLVALUE) {
      meanY = C_adj * xtt1 + DwReg_adj;
      VarY  = C_adj * Ptt1 * C_adj.t() + R_adj;
      //
      // 
      //
      //
      //
      //
      // chol_VarY = arma::chol(VarY, "lower");
      // dim_tmp = chol_VarY.n_cols;
      // chol_VarY_diag = chol_VarY.diag();
      // arma::vec vec_tmp_ll(dim_tmp);
      // 
      // for (unsigned int dd = 0; dd < dim_tmp; ++dd) {
      //   // for(dd = 0; dd < d; dd++) {
      //   acc = 0.0;
      //   for (unsigned int ii = 0; ii < dd; ii++) acc += vec_tmp_ll(ii) * chol_VarY.at(dd, ii);
      //   vec_tmp_ll(dd) = (yObs_adj(dd) - meanY(dd) - acc) / chol_VarY_diag(dd);
      // }
      // ll_out_tmp += arma::sum(arma::square(vec_tmp_ll));
      //
      // 
      //
      //
      //
      //
      ll_out_tmp += log(arma::det(VarY));
      ll_out_tmp += arma::sum((yObsAdj - meanY).t() * arma::inv(VarY) * (yObsAdj - meanY));
    }
  // 
  // part1 <- -(TT*dimY - MissingsObsTotal) /2 * log(2*pi)
  //   
  //   llOUT <- part1 - 0.5 * part2
  // if(isFALSE(LOG)) {
  //   llOUT <- exp(llOUT)
  // } else if(!isTRUE(LOG) || isFALSE(LOG)) {
  //   stop("Wrong argument type for 'LOG': must be logical either TRUE or FALSE.")
  // }
  // return(llOUT)
  // period t+1 quantities for next iteration
    if (PDSTORE) {
      xtt1STORE.col(t + 1)   = xtt1;
      Ptt1STORE.slice(t + 1) = Ptt1;
    }
  }
  if (PDSTORE & !LLVALUE) {
    return(Rcpp::List::create(Rcpp::Named("mfdEXP") = xtt,
                              Rcpp::Named("mfdVAR") = Ptt,
                              Rcpp::Named("pddEXP") = xtt1STORE,
                              Rcpp::Named("pddVAR") = Ptt1STORE));
    
  } if (LLVALUE) {
    // ll_out = ll_out_tmp;
    // ll_out = -0.5 * ll_out_tmp - 0.5 * (TT * dimY - missing_obs_total) * log(2.0 * arma::datum::pi);
    //
    //
    //
    //
    //
    ll_out = -(TT*dimY - missing_obs_total) /2 * log(2.0 * arma::datum::pi) - 0.5 * ll_out_tmp;
    return(Rcpp::List::create(Rcpp::Named("mfdEXP") = xtt,
                              Rcpp::Named("mfdVAR") = Ptt,
                              Rcpp::Named("kfLLH") = ll_out));
  } else {
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
arma::mat bs(const arma::uword TT,
             const int nfac,
             const arma::mat& Phi,
             const arma::mat& Q,
             const arma::mat& filt_f,
             const arma::cube& filt_P) {
  int TT_use = static_cast<int>(TT);
  arma::mat fPost(nfac, TT_use);

  arma::colvec fTT = filt_f.col(TT_use - 1);
  arma::colvec ftt(TT_use);

  arma::mat PTT = filt_P.slice(TT_use - 1);
  arma::mat Ptt(arma::size(PTT));

  fPost.col(TT_use - 1) = arma::mvnrnd(fTT, PTT);
  
  arma::mat IS(TT_use, TT_use);
  arma::mat sim_v(TT_use, TT_use);
  arma::colvec sim_m(TT_use);

  for (int t = TT_use - 2;  t >= 0; --t) {
    ftt = filt_f.col(t);
    Ptt = filt_P.slice(t);
    IS = Ptt * Phi.t() * arma::inv(Phi * Ptt * Phi.t() + Q);
    sim_m = ftt + IS * (fPost.col(t + 1) - Phi * ftt);
    sim_v = Ptt - IS * Phi * Ptt;
    sim_v = 0.5 * (sim_v + sim_v.t());

    fPost.col(t) = arma::mvnrnd(sim_m, sim_v);
  }
  return(fPost);
}