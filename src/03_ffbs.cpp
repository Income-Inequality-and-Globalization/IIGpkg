#include "03_ffbs.h"
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
                 bool PDSTORE) {
  // xtt <- matrix(0, nrow = dimX, ncol = TT)
  arma::mat xtt(dimX, TT);
  // Ptt <- array(0, dim = c(dimX, dimX, TT),
               // dimnames = list(NULL, NULL, as.character(1:TT)))
  arma::cube Ptt(dimX, dimX, TT);
  if (PDSTORE) {
    // xtt1STORE <- matrix(0, nrow = dimX, ncol = TT + 1)
    arma::mat xtt1STORE(dimX, TT + 1);
    // Ptt1STORE <- array(0, dim = c(dimX, dimX, TT + 1),
    //                    dimnames = list(NULL, NULL, as.character(1:(TT + 1))))
    arma::cube Ptt1STORE(dimX, dimX, TT + 1);
  }
  // BuRegInit <- computeMatReg(mat = B, reg = u00, dim = dimX, lenT = 1)
  // BuReg     <- computeMatReg(mat = B, reg = uReg, dim = dimX, lenT = TT)
  // DwReg     <- computeMatReg(mat = D, reg = wReg, dim = dimY, lenT = TT)
  arma::mat D_wreg = compute_mat_reg(D, wReg);
// 
// 
//     xtt1 <- computeXtt1(A, x00, BuRegInit[, 1], dimX)
//     Ptt1 <- computePtt1(A, P00, Q)
// 
//     if (isTRUE(PDSTORE)) {
//       xtt1STORE[, 1]   <- xtt1
//       Ptt1STORE[, , 1] <- Ptt1
//     }
// 
// 
//     for(t in 1:TT) {
//       AllObsMissing <- FALSE
//       CAdj <- C
//       if(is.na(dim(R)[3])){
//         RAdj <- R
//       }else{
//         RAdj <- R[, , t]
//       }
//       yObsAdj <- yObs[, t]
//       DwRegAdj <- DwReg[, t]
//       MissingObs <- which(is.na(yObsAdj))
// 
//         if (length(MissingObs) != 0) {
//           if (length(MissingObs) == dimY) {
//             AllObsMissing <- TRUE
//           } else {
//             W <- diag(dimY)
//             W <- W[-MissingObs, ]
// 
//             yObsAdj[MissingObs] <- 0
//             CAdj[MissingObs, ] <- 0
//             RAdj[MissingObs, ] <- 0
//             RAdj[ ,MissingObs] <- 0
// 
//             CAdj <- W %*% CAdj
//             RAdj <- W %*% RAdj %*% t(W)
//             yObsAdj <- W %*% yObsAdj
//             DwRegAdj <- W %*% DwRegAdj
//           }
//         }
// 
//         if(!AllObsMissing) {
//           Lt   <- computeLt(CAdj, Ptt1, RAdj)
//           Kt   <- computeKt(Ptt1, CAdj)
// 
// 
// # period t quantities for current iteration
//           kGain      <- computekG(yObsAdj, CAdj, xtt1, DwRegAdj)
//             xtt[, t]   <- computeXtt(xtt1, Kt, Lt, kGain)
//             Ptt[, , t] <- computePtt(Ptt1, Kt, Lt, CAdj, RAdj, t)
//         } else {
//           xtt[, t] <- xtt1
//           Ptt[, , t] <- Ptt1
//         }
// 
// 
//         xtt1 <- computeXtt1(A, xtt[, t], BuReg[, t], dimX)
//           Ptt1 <- computePtt1(A, Ptt[, , t], Q)
// 
// # period t+1 quantities for next iteration
//           if (isTRUE(PDSTORE)) {
//             xtt1STORE[, t + 1]   <- xtt1
//             Ptt1STORE[, , t + 1] <- Ptt1
//           }
// 
//     }
//     if (PDSTORE) {
//       out <- list(mfdEXP = xtt, mfdVAR = Ptt,
//                   pddEXP = xtt1STORE, pddVAR = Ptt1STORE)
// 
//     } else {
//       out <- list(mfdEXP = xtt, mfdVAR = Ptt)
//     }
//     return(out)
  return(Rcpp::List::create(Rcpp::Named("a") = 1));
}