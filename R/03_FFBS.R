compute_FFBS <- function(yObs, uReg, wReg,
                         num_fac, num_y, N_num_y, 
                         TT, NN,
                         VhatArray_A, Phi,
                         B, C, D, Q, 
                         initX, initU, initP, 
                         PDSTORE, try_catch_errors = NULL) {
  # if (!ident_block) {
  # Vhat Array (mit Adjustmentmatrix A) bzgl. der Zeit sortiert
  R <- bdiagByTime(
    VhatArray_A = VhatArray_A,
    npara = num_y,
    N = NN,
    TT = TT,
    Nnpara = N_num_y
  )
  if (is.null(try_catch_errors)) {
    # Kalman-Filer
    KF <- RcppSMCkalman::kfMFPD(
      yObs = yObs, uReg = uReg, wReg = wReg,
      dimX = num_fac, dimY = N_num_y, TT = TT,
      A = Phi, B = B, C = C, D = D,
      Q = Q, R = R, x00 = initX,
      u00 = initU, P00 = initP, PDSTORE = PDSTORE)
    # Backward Sampling and return
    return(GibbsSSM_f(TT = TT, nfac = num_fac,
                      Phi = Phi, Q = Q,
                      filt_f = KF$mfdEXP, filt_P = KF$mfdVAR))
  } else {
    stopifnot(`Arg. 'try_catch_errors' must be a named list` = 
      all(names(try_catch_errors) %in% c("fSTORE", " BSTORE", "ASTORE",
                                         "block_count", "initials",
                                         "storePath_adj", "storePath_kfe",
                                         "store_count", "iter")))
    invisible(
      capture.output(
        KF <- tryCatch(
          {
            RcppSMCkalman::kfMFPD(
              yObs = yObs, uReg = uReg, wReg = wReg,
              dimX = num_fac, dimY = N_num_y, TT = TT,
              A = Phi, B = B, C = C, D = D,
              Q = Q, R = R, x00 = initX,
              u00 = initU, P00 = initP, PDSTORE = PDSTORE
            )
          },
          error = function(e) {e$message}
        )
      )
    )
    # Speichert moegliche Fehlermeldung des Kalman-Filter und die Zwischenergebnisse zum Zeitpunkt des Fehlers zurueck
    if (is.character(KF)) {
      msg_error_kf <- KF
      if (try_catch_errors$store_count == 0) {
        dir.create(try_catch_errors$storePath_adj, recursive = TRUE)
      }
      storePath_kf <- get_store_path_kf_error_iter(
        try_catch_errors$storePath_kfe, try_catch_errors$iter
      )
      out_error_kf <- list(
        f = try_catch_errors$fSTORE,
        B = try_catch_errors$BSTORE)
      if (isTRUE(VdiagEst)) {
        out_error_kf <- c(out_error_kf, list(V = try_catch_errors$VSTORE))
      } else if (isFALSE(VdiagEst)) {
        out_error_kf <- c(out_error_kf, list(A = try_catch_errors$ASTORE))
      }
      out_error_kf <- c(out_error_kf,
                        list(
                          blockCount = try_catch_errors$block_count,
                          errorMsg = msg_error_kf,
                          initials = try_catch_errors$initials)
                        )
      saveRDS(out_error_kf, file = storePath_kf)
      stop("KF produced numerical errors.")
    }
    return(GibbsSSM_f(TT = TT, nfac = num_fac, Phi = Phi, Q = Q,
                      filt_f = KF$mfdEXP, filt_P = KF$mfdVAR))
    # return(list(f = try_catch_errors$fSTORE,
    # B = try_catch_errors$BSTORE,
    # D = try_catch_errors$DSTORE,
    # A = ASTORE, 
    # blockCount = block_count, errorMsg = msg_error_kf))
    #
    #
    #
    #
    #
    # ident_block <- FALSE # wird nicht mehr benoetigt, nicht auskommentiert, da es unten auch noch drin steht. War vom alten Sampler, um die identifizierenden Restriktionen einzuhalten.
    # Backward Sampling
    # try({KF <-RcppSMCkalman::kfMFPD(yObs = yObs, uReg = uReg, wReg = wReg,
    #                                                      dimX = num_fac, dimY = Nnpara, TT = TT,
    #                                                      A = Phi, B = NULL, C = B, D = c_,
    #                                                      Q = Q, R = VhatArrayBdiagByTime, x00 = initX,
    #                                                     u00 = initU, P00 = initP, PDSTORE = F)})
    # if(inherits(KF, "try-error")){
    #   print("Fehler")
    #   break
    # }
    # }
  }
}
get_store_path_kf_error_base <- function(VdiagEst, store_path_adj, init_set, inc_obs_new) {
  if (VdiagEst) {
  paste0(
    store_path_adj,
    "/",
    "B",
    round(init_set$B0[1, 1, 1], 2),
    "_Omega", init_set$Omega0[1, 1],
    "_V", init_set$Vhat[1, 1, 1],
    "_alpha", init_set$alpha0,
    "_beta", init_set$beta0,
    "_IO", inc_obs_new, "_error"
   )
  } else {
    paste0(
      store_path_adj, "/",
      "B", round(init_set$B0[1, 1, 1], 2),
      "_Omega", init_set$Omega0[1, 1],
      "_A", init_set$A[1, 1],
      "_Psi", init_set$Psi0[1, 1],
      "_nu0", init_set$nu0,
      "_IO", inc_obs_new, "_error"
    )
  }
}
get_store_path_kf_error_iter <- function(store_path_kf_error_base, iter) {
  paste0(store_path_kf_error_base, "_ITER_",iter, ".rds")
}
#' Backward Sampling based Kalman-Filter Output
#'
#' @param TT time
#' @param nfac Number of factors
#' @param Phi phi param
#' @param Q VCM part
#' @param filt_f Means of filtering distribution from Kalman-Filter
#' @param filt_P Variances of filtering distribution from Kalman-Filter
#'
#' @return posterior filtered states
#' @export
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
