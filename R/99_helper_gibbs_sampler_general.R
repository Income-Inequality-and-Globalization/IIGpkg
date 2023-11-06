set_simulate_factors <- function(fPost) {
  if (missing(fPost)) {
    return(TRUE)
  }
  FALSE
}
set_reg_specification <- function(wRegSpec) {
  if (missing(wRegSpec)) wRegSpec <- 0
  wRegSpec
}
set_store_path_subdir <- function(store_path, V_DIAG_EST, SAMPLE_A,
                                  init_set, V_start, cov_scale, num_joint_fac,
                                  w_reg_spec, Omega_D0, inc_obs_new) {
  if (store_path == "none") {
    return(invisible(list(store_path_adj = "none",
                          store_path_rds = "none")))
  }
  base_adj <- paste0(
    store_path,
    "/", "cs",
    cov_scale,
    "_pj", num_joint_fac,
    "_Reg", w_reg_spec,
    "_B", round(init_set$B0[1, 1, 1], 2),
    "_Om",
    init_set$Omega0[1, 1], "_D",
    init_set$D0[1, 1, 1], "_OmD",
    Omega_D0[dim(Omega_D0)[1], dim(Omega_D0)[1]]
  )
  # Zu `Omega_D0[dim(Omega_D0)[1], dim(Omega_D0)[1]]`:
  # -> speichert Startwert der Priorvarianz fuer die
  # Makroregressor-Koeffizienten. Da ich Omega0 immer diagonal waehle und fuer
  # jeden Makroregressor-Koeffizienten die gleiche Varianz waehle, reicht es
  # hier einen Parameter zu speichern.
  if (store_path != "none" & !missing(store_path)) {
    if (V_DIAG_EST || SAMPLE_A) {
      store_path_adj <- paste0(
        base_adj,
        "_A", init_set$A[1, 1],
        "_Psi", init_set$Psi0[1, 1],
        "_nu", init_set$nu0,
        "_IO", inc_obs_new
      )
    } else {
      store_path_adj <- paste0(
        base_adj,
        "_IO",
        inc_obs_new
      )
    }
  }
  # if (is.null(store_path_adj)) stop("HIT ME IN MY FACE")
  store_path_rds <- get_store_path_rds(store_path_adj,
                                       V_DIAG_EST,
                                       SAMPLE_A,
                                       init_set,
                                       V_start,
                                       cov_scale,
                                       num_joint_fac,
                                       Omega_D0,
                                       inc_obs_new)
  store_path_omg <- paste0(
    store_path_adj, 
    "/", "invOmega1_", "pj", num_joint_fac, 
    "_B", round(init_set$B0[1, 1, 1], 2),
    "_Omega", init_set$Omega0[1, 1],
    "_A", init_set$A[1, 1],
    "_Psi", init_set$Psi0[1, 1],
    "_nu0", init_set$nu0,
    "_IO", inc_obs_new)
  # STORE WITH uSTORE:
  # saveRDS(list(f = fSTORE, B = BSTORE, D = DSTORE, V = VSTORE, u = uSTORE, blockCount = block_count,  errorMsg = errorMsg, initials = initials)
  #         , file = paste0(storePath_adj,"/", "pj",njointfac,"_B",round(initials$B0[1,1,1],2),"_Om",initials$Omega0[1,1],"_D",initials$D[1,1,1],"_OmD",Omega0[dim(Omega0)[1], dim(Omega0)[1]],"_V",Vstart, "_alpha",initials$alpha0,"_beta",initials$beta0,"_IO",incObsNew,".rds") )
  return(list(store_path_adj = store_path_adj,
              store_path_rds = store_path_rds,
              store_path_omg = store_path_omg))
}
get_store_path_rds <- function(store_path_adj,
                               V_DIAG_EST, SAMPLE_A,
                               init_set, V_start,
                               cov_scale, num_joint_fac, Omega_D0,
                               inc_obs_new) {
  base_path <- file.path(store_path_adj,
                         paste0("cs", cov_scale,
                                "_pj", num_joint_fac,
                                "_B", round(init_set$B0[1, 1, 1], 2),
                                "_Om", init_set$Omega_D0[1, 1],
                                "_D", init_set$D0[1, 1, 1],
                                "_OmD", Omega_D0[dim(Omega_D0)[1], dim(Omega_D0)[1]])
  )
  if (V_DIAG_EST) {
    base_path <- paste0(base_path,
                        "_V", V_start,
                        "_alpha", init_set$alpha0,
                        "_beta", init_set$beta0)
  } else if (SAMPLE_A) {
    base_path <- paste0(base_path,
                       "_A", init_set$A[1, 1],
                       "_Psi", init_set$Psi0[1, 1],
                       "_nu", init_set$nu0)
  }
  paste0(base_path, "_IO", inc_obs_new, ".rds")
}
store_mcmc <- function(VdiagEst, init_set,
                       fSTORE, BSTORE, DSTORE, VSTORE, ASTORE,
                       block_count, msg_error_kf, store_path_rds) {
  if (store_path_rds == "none") return(invisible(NULL))
  store_mcmc_rds <- list(f = fSTORE,
                         B = BSTORE,
                         D = DSTORE)
  if (VdiagEst) {
    store_mcmc_rds <- c(store_mcmc_rds, list(V = VSTORE))
  } else {
    store_mcmc_rds <- c(store_mcmc_rds, list(A = ASTORE))
  }
  saveRDS(c(store_mcmc_rds,
            list(blockCount = block_count,
                 errorMsg = msg_error_kf,
                 initials = init_set)),
          file = store_path_rds)
}
get_check_store <- function(store_path, mcmc_iter, mcmc_itermax, store_unit) {
  if (store_path == "none") return(FALSE)
  CHECK_ITER_REACH <- (mcmc_iter == 100)
  CHECK_ITER_MAX   <- (mcmc_iter == mcmc_itermax)
  CHECK_STORE_UNIT <- (mcmc_iter %% store_unit == 0)
  return(CHECK_ITER_REACH || CHECK_ITER_MAX || CHECK_STORE_UNIT)
}
get_fails <- function(out_gibbs, iter_fail) {
  iter_fail_before <- iter_fail - 1
  
  f_prior_to_fai <- out_gibbs$Gibbs2_SM_SA$f[, , iter_fail_before]
  B_prior_to_fai <- out_gibbs$Gibbs2_SM_SA$B[, , iter_fail_before]
  D_prior_to_fai <- out_gibbs$Gibbs2_SM_SA$D[, , iter_fail_before]
  A_prior_to_fai <- out_gibbs$Gibbs2_SM_SA$A[, , iter_fail_before]
  V_prior_to_fai <- out_gibbs$Gibbs2_SM_SA$V[, iter_fail_before]
  
  return(list(f_prior_to_fai = f_prior_to_fai,
              B_prior_to_fai = B_prior_to_fai,
              D_prior_to_fai = D_prior_to_fai,
              A_prior_to_fai = A_prior_to_fai,
              V_prior_to_fai = V_prior_to_fai))
}