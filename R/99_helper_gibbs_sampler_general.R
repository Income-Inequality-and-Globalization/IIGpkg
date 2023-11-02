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
                                  init_set, cov_scale, num_joint_fac,
                                  w_reg_spec, Omega_D0, inc_obs_new) {
  if (store_path == "none") {
    return(invisible(store_path))
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
  if (is.null(store_path_adj)) stop("HIT ME IN MY FACE")
  return(store_path_adj)
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