#include "03_ffbs_helper.h"
arma::mat compute_mat_reg(const arma::mat& mat, 
                          const arma::mat& reg) {
  // int dim,
  // int lenT) {
  // if (is.null(mat) || is.null(reg)) return(matrix(0, nrow = dim, ncol = lenT))
  return(mat * reg);
}
arma::mat copute_Xtt_1(const arma::mat& A, const arma::mat& xtt) {
  // , BuReg, dimX) {
  // if (dimX == 1) return(A*xtt + BuReg)
  return(A * xtt);
  // computeXtt1 <- function(A, xtt, BuReg, dimX) {
  //   if (dimX == 1) return(A*xtt + BuReg)
  //     A %*% xtt + BuReg
  // }
}
arma::mat compute_ptt_1(const arma::mat& A,
                        const arma::mat&Ptt,
                        const arma::mat&Q) {
  // computePtt1 <- function(A, Ptt, Q) {
  //   A %*% tcrossprod(Ptt, A) + Q
  //   A %*% P %*% t(A) + Q
  // }
  return(A * Ptt * A.t() + Q);
}
arma::mat compute_Lt(const arma::mat& C,
                     const arma::mat& Ptt1,
                     const arma::mat& R) {
  // computeLt <- function(C, Ptt1, R) {
  // #solve(C %*% tcrossprod(Ptt1, C) + R)
  //   R <- 0.5*(R + t(R))
  //   R_inv <- solve(R)
  //   P <- Ptt1
  //   P <- 0.5*(P + t(P))
  //   P_inv <- solve(P)
  //   PCRC <- solve(P_inv + t(C) %*% R_inv %*% C)
  //   R_inv - R_inv %*% C %*% PCRC %*% t(C) %*% R_inv
  // }
  arma::mat R2 = 0.5 * (R + R.t());
  arma::mat R_inv = arma::inv(R);
  arma::mat P = 0.5 * (Ptt1 + Ptt1.t());
  arma::mat P_inv = arma::inv(P);
  
  arma::mat PCRC = arma::inv(P_inv + C.t() * R_inv * C);
  return(R_inv - R_inv * C * PCRC * C.t() * R_inv);
}
arma::mat compute_Kt(const arma::mat& Ptt1, const arma::mat& C) {
  // computeKt <-function(Ptt1, C) {
  //   tcrossprod(Ptt1, C)
  // }
  return(Ptt1 * C.t());
}
arma::mat compute_kG(const arma::mat& yObs, const arma::mat& C,
                     const arma::mat& xtt1, const arma::mat& DwReg) {
  // computekG <- function(yObs, C, xtt1, DwReg) {
  //   yObs - C %*% xtt1 - DwReg
  // }
  return(yObs - C * xtt1 - DwReg);
}
arma::mat compute_Xtt(const arma::mat& xtt1, const arma::mat& Kt,
                      const arma::mat& Lt, const arma::mat& kGain) {
  // computeXtt <- function(xtt1, Kt, Lt, kGain) {
  //   xtt1 + Kt %*% Lt %*% kGain
  // }
  return(xtt1 + Kt * Lt * kGain);
}
arma::mat compute_Ptt(const arma::mat& Ptt1, const arma::mat& Kt,
                      const arma::mat& Lt, const arma::mat& C, 
                      const arma::mat& R) {
  // computePtt <- function(Ptt1, Kt, Lt, C, R, t) {
  // #Ptt <- Ptt1 - Kt %*% tcrossprod(Lt, Kt)
  //   I <- diag(dim(Ptt1)[1])
  //   KL <- Kt %*% Lt
  //   IKLC <- I - KL %*% C
  //   Ptt <- IKLC %*% tcrossprod(Ptt1, IKLC) + KL %*% tcrossprod(R, KL)
  // 
  // # if (!matrixcalc::is.positive.definite(Ptt)) {
  // # stop(paste0("matrix is no longer p.d. at iteration number: ", t))
  // # }
  //   return(Ptt)
  // }
  arma::mat I = arma::eye(arma::size(Ptt1));
  arma::mat KL = Kt * Lt;
  arma::mat IKLC = I - KL * C;
  return(IKLC * Ptt1 * IKLC.t() + KL * R * KL.t());
}