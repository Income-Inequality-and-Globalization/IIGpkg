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
set_store_path_subdir <- function(store_path, SAMPLE_V, SAMPLE_A, SAMPLE_H,
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
    store_path_adj <- base_adj
    if (SAMPLE_V || SAMPLE_A) {
      store_path_adj <- paste0(
        store_path_adj,
        "_A", init_set$A[1, 1],
        "_Psi", init_set$Psi0[1, 1],
        "_nu", init_set$nu0
      )
    } 
    if (SAMPLE_H) {
      store_path_adj <- file.path(
        store_path_adj, 
        paste0(
          "mu_b0_", init_set$prior_list$hyperpriors$mu_b0[1], "_",
          "Sigma_b0_", init_set$prior_list$hyperpriors$Sigma_b0[1, 1], "_",
          "alpha_b0_", init_set$prior_list$hyperpriors$alpha_b0[1], "_",
          "beta_b0_", init_set$prior_list$hyperpriors$beta_b0[1]
        )
      )
    }
    store_path_adj <- paste0(
      store_path_adj,
      "_IO",
      inc_obs_new
    )
  }
  # if (is.null(store_path_adj)) stop("HIT ME IN MY FACE")
  store_path_rds <- get_store_path_rds(store_path_adj,
                                       SAMPLE_V,
                                       SAMPLE_A,
                                       SAMPLE_H,
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
  store_path_kfe <- get_store_path_kf_error_base(SAMPLE_V,
                                                 store_path_adj,
                                                 init_set, 
                                                 inc_obs_new)
  return(list(store_path_adj = store_path_adj,
              store_path_rds = store_path_rds,
              store_path_omg = store_path_omg,
              store_path_kfe = store_path_kfe)
         )
}
get_store_path_rds <- function(store_path_adj,
                               SAMPLE_V, SAMPLE_A, SAMPLE_H,
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
  if (SAMPLE_V) {
    base_path <- paste0(base_path,
                        "_V", V_start,
                        "_alpha", init_set$alpha0,
                        "_beta", init_set$beta0)
  }
  if (SAMPLE_A) {
    base_path <- paste0(base_path,
                       "_A", init_set$A[1, 1],
                       "_Psi", init_set$Psi0[1, 1],
                       "_nu", init_set$nu0)
  }
  if (SAMPLE_H) {
    base_path <- paste0(
      base_path,
      "mu_b0_", init_set$prior_list$hyperpriors$mu_b0[1], "_",
      "Sigma_b0_", init_set$prior_list$hyperpriors$Sigma_b0[1, 1], "_",
      "alpha_b0_", init_set$prior_list$hyperpriors$alpha_b0[1], "_",
      "beta_b0_", init_set$prior_list$hyperpriors$beta_b0[1]
    )
  }
  paste0(base_path, "_IO", inc_obs_new, ".rds")
}
store_mcmc <- function(output_list, store_path_rds) {
  if (store_path_rds == "none") return(invisible(NULL))
  saveRDS(output_list, file = store_path_rds)
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
get_initials <- function(envir_list) {
  init_names <- c("itermax",
    "identmax",
    "npara",
    "nreg",
    "njointfac",
    "yObs",
    "c_",
    "D",
    "D0",
    "B",
    "Phi",
    "Q",
    "type",
    "initX",
    "initP",
    "initU",
    "wRegSpec",
    "wReg",
    "uReg",
    "B0",
    "Omega0",
    "selectR",
    "selectC",
    "Vhat",
    "incObsOld",
    "incObsNew",
    "covScale",
    "VhatDiagScale",
    "VhatDiagScale_start",
    "VdiagEst",
    "alpha0",
    "beta0",
    "countryA",
    "A",
    "scaleA",
    "diagA",
    "Psi0",
    "nu0",
    "shape0",
    "rate0",
    "storePath",
    "fPost",
    "sampleA",
    "identification",
    "storeUnit"
  )

  prior_list <- envir_list[["prior_list"]]
  remove_names <- c("prior_list")
  envir_list[remove_names] <- NULL
  envir_list_out <- envir_list[1:19]
  envir_list_out$B0     <- prior_list$priors$B0
  envir_list_out$Omega0 <- prior_list$priors$Omega0
  envir_list_out <- c(envir_list_out, envir_list[20:43])
  # B0       <- prior_list$priors$B0
  # mu_b0    <- prior_list$hyperpriors$mu_b0_par
  # Sigma_b0 <- prior_list$hyperpriors$Sigma_b0
  # Omega0   <- prior_list$priors$Omega0
  # alpha_b0 <- prior_list$hyperpriors$alpha_b0
  # beta_b0  <- prior_list$hyperpriors$beta_b0
  return(envir_list_out)
}
plot_mcmc_ffbs <- function(f_stored, mm, 
                           f_trues = NULL,
                           options_f =list(
                             grid = c(3,2),
                             id_f = 1:6),
                           options_p = list(
                             break_mm = 50,
                             pause = 0.4)
                           ) {
  if (mm != 1 && (mm %% options_p$break_mm != 0)) return(invisible(NULL))
  f_test <- f_stored[, ,  ceiling(mm / 2):mm]
  f_test_mean <- apply(f_test, c(1,2), mean)
  f_test_sd   <- apply(f_test, c(1,2), sd)
  if (is.na(f_test_sd[1,1])) f_test_sd[] <- 1

  par(mfrow = options_f$grid)
  for (tt in options_f$id_f) {
    t_taken_loc <- tt
    ci_lower <- f_test_mean[t_taken_loc, ] - 1.96*f_test_sd[t_taken_loc, ]
    ci_upper <- f_test_mean[t_taken_loc, ] + 1.96*f_test_sd[t_taken_loc, ]
    if (!is.null(f_trues)) {
      ylim_min <- min(min(ci_lower), min(f_trues[t_taken_loc, ]))
      ylim_max <- max(max(ci_upper), max(f_trues[t_taken_loc, ]))
    } else {
      ylim_min <- min(ci_lower)
      ylim_max <- max(ci_upper)
    }
    plot(f_test_mean[t_taken_loc, ], type = "l",
         ylim = c(ylim_min, ylim_max),
         ylab = paste0("fac. No. ", tt), xlab = "t",
         col = "blue")
    if (!is.null(f_trues)) {
      lines(f_trues[t_taken_loc, ], col = "green")
    }
    lines(ci_lower, col = "blue", lty = "dashed")
    lines(ci_upper, col = "blue", lty = "dashed")
  }
  mtext("Factors ----- blue: filtered wit 95%-CI ----- green: true",
        side = 3, line = -2, outer = TRUE)
  par(mfrow = c(1, 1))
  if (!is.null(options_p$pause)) Sys.sleep(options_p$pause)
  return(invisible(NULL))
}
