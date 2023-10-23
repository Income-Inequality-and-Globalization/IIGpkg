#' Backward Sampling based Kalman-Filter Output
#' 
#' @inheritParams GibbsSSM_2
#' @param TT
#' @param nfac Number of factors
#' @param Phi
#' @param Q
#' @param filt_f Means of filtering distribution from Kalman-Filter 
#' @param filt_P Variances of filtering distribution from Kalman-Filter
#'
#' @return
#' @export
#'
#' @examples
GibbsSSM_f <- function(TT, nfac, Phi, Q, filt_f, filt_P) {
  fTT <- filt_f[, TT]
  PTT <- filt_P[, , TT]

  fPost <- matrix(rep(0, nfac * TT), ncol = TT)
  fPost[, TT] <- MASS::mvrnorm(n = 1, mu = fTT, Sigma = PTT)

  for (t in (TT - 1):1) {
    ftt <- filt_f[, t]
    Ptt <- filt_P[, , t]
    IS <- Ptt %*% t(Phi) %*% solve(Phi %*% Ptt %*% t(Phi) + Q)
    sim_m <- ftt + IS %*% (fPost[, t + 1] - Phi %*% ftt)
    sim_v <- Ptt - IS %*% Phi %*% Ptt
    sim_v <- 0.5 * (sim_v + t(sim_v))
    # if(!matrixcalc::is.positive.definite(sim_v)){
    #   print(iter)
    #   print(sim_v)
    #   print(filt_P)
    # }
    fPost[, t] <- MASS::mvrnorm(n = 1, mu = sim_m, Sigma = sim_v)
  }

  return(fPost)
}
