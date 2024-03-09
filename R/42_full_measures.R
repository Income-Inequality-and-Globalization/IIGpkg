#' Precomputes large container that allow to get-posterior measures
#'
#' @param out_gibbs the Gibbs output as an object of class `GibbsOutputIIG`
#' @param transformation_infos the centering and standardizing vectors for the
#'    `y-measurement` data as a vector for all countries (long index) and all
#'    parameters (short index) e.g. a vector of length `30` (`3 x 10`, for ten
#'    countries)
#' @param const_num_para numeric value giving the number of parameters in the
#'    model
#'
#' @return a named list of precomputed container and regressors/regressor grids
#'   used to get them
#' @export
precompute_measures_me <- function(out_gibbs,
                                   transformation_infos,
                                   const_num_para = 3) {
  if (missing(transformation_infos)) stop("Arg. missing.")
  validate_GibbsOutputIIG(out_gibbs)

  names_regs  <- transformation_infos$x_settings$names_regs
  cut_off_num <- transformation_infos$x_settings$cut_off_num

  NN_num_par <- dim(out_gibbs$B)[1]
  NN <- NN_num_par / const_num_para
  TT <- dim(out_gibbs$f)[2]
  MM <- dim(out_gibbs$B)[3]
  KK <- length(names_regs)
  GG <- transformation_infos$x_settings$grid_length

  B   <- get_cnt_me_B(out_gibbs$B)
  f   <- get_cnt_me_f(out_gibbs$f)
  D   <- get_cnt_me_D(out_gibbs$D, names_regs)
  # D_u <- get_cnt_me_D_uncertainty(D, names_regs, KK, MM)

  wRegs     <- get_cnt_wRegs(transformation_infos$meta$wRegs_grd, names_regs)
  wRegs_raw <- get_cnt_wRegs(transformation_infos$meta$wRegs_raw, names_regs)
  wRegs_plus_sd <- get_cnt_fd_wRegs(transformation_infos$meta$wRegs_grd,
                                    names_regs)

  wRegs_grid         <- get_cnt_WR(wRegs, names_regs, GG, cut_off_num)
  wRegs_grid_raw     <- get_cnt_WR(wRegs_raw, names_regs, GG, cut_off_num)

  # B_m <- apply(B, c(1, 2), mean)
  # f_m <- apply(f, c(1, 2), mean)

  post_par_estimates_fd <- get_cnt_par_marg_estim(NN = NN, TT = TT, MM = MM, names_regs)
  post_par_estimates    <- get_cnt_par_post_estim(NN = NN, TT = TT, MM = MM)
  post_me_grid          <- get_cnt_par_meas_estim(NN, TT, GG, MM, names_regs)
  # post_me_raws   <- get_cnt_par_me_raws(NN, TT, MM, names_regs)
  for (mm in seq_len(MM)) {
    Bxf <- B[, , mm] %*% f[, , mm]
    post_par_estimates[, , mm] <- Bxf + D[, , mm] %*% wRegs
    for (kk in seq_len(KK)) {
      # post_me_raws[, , mm, kk] <- Bxf + D_u[, , mm, kk] %*% wRegs
      post_par_estimates_fd[, , mm, kk] <- Bxf + D[, , mm] %*% wRegs_plus_sd[, , kk]
      for (gg in seq_len(GG)) {
        # post_me_grid[, , gg,  mm, kk] <- Bxf + D_u[, , mm, kk] %*% wRegs_grid[, , gg, kk]
        post_me_grid[, , gg,  mm, kk] <- Bxf + D[, , mm]  %*% wRegs_grid[, , gg, kk]
        # post_me[, , gg,  mm, kk] <- (B_m %*% f_m +
        #  D_unc[ , , mm, kk] %*% WR[, , gg, kk])
      }
    }
    progress_any(mm, MM)
  }
  return(list(post_par_estimates = post_par_estimates,
              post_par_estimates_forward_diff = post_par_estimates_fd,
              post_me_grid = post_me_grid,
              # post_me_raws = post_me_raws,
              regressor_grid = list(
                reg_grid_transformed = wRegs_grid,
                reg_grid_original = wRegs_grid_raw)
              )
  )
}
#' Generates various posterior measures from precomputed results
#'
#' Precomputed results are obtained from [precompute_measures_me]
#'
#' @param par_post_estim estimated posterior fitted values container
#' @param post_me estimated posterior fitted values container additionally
#'    generated on a grid and per regressor
#' @param regs_to_use a vector of regressors to compute ME measures for
#' @inheritParams precompute_measures_me
#'
#' @return a named list of measures
#'
#' @export
generate_measures_me <- function(par_post_estim,
                                 par_post_estim_fd,
                                 post_me,
                                 regs_to_use,
                                 transformation_infos,
                                 const_num_para = 3) {
  stopifnot(dim(post_me) != 5)
  dim_last_meas <- 5
  dim_last_marg <- 4
  dim_post <- 3
  dim_marg <- 3
  dim_meas <- 4

  scale_vals  <- transformation_infos$y_settings$scaling
  centr_vals  <- transformation_infos$y_settings$centering
  names_regs  <- transformation_infos$x_settings$names_regs


  par_post_estim_bt    <- f_bt_pars(par_post_estim, centr_vals, scale_vals)

  a_post_bt <- get_subset_par(par_post_estim_bt, par_name = "a", dim_post)
  q_post_bt <- get_subset_par(par_post_estim_bt, par_name = "q", dim_post)
  m_post_bt <- get_subset_par(par_post_estim_bt, par_name = "mu", dim_post)

  mu_info_TT      <- compute_mu_info(m_post_bt)
  gini_info_TT    <- compute_gini_info(a_post_bt, q_post_bt)

  # mu_info_TT_2   <- compute_mu_info(m_post_bt, 0.8)
  # gini_info_TT_2 <- compute_gini_info(a_post_bt, q_post_bt, 0.8)

  KK <- length(regs_to_use)
  mu_info_KK      <- get_cnt_info_KK(KK, regs_to_use)
  gini_info_KK    <- get_cnt_info_KK(KK, regs_to_use)
  mu_info_fd_TT   <- get_cnt_info_KK(KK, regs_to_use)
  gini_info_fd_TT <- get_cnt_info_KK(KK, regs_to_use)

  regs_to_use_id <- which(names_regs %in% regs_to_use)
  post_me <- post_me[, , , , regs_to_use_id, drop = FALSE]
  for (kk in 1:KK) {
    cat(paste0("Iteration of KK-type measures -- kk = ", kk, "\n"))
    post_meas_kk <- abind::adrop(post_me[, , , , kk, drop = FALSE], dim_last_meas)
    tmp_meas_bt  <- f_bt_meas(post_meas_kk, centr_vals, scale_vals)

    a_post_me_KK_bt <- get_subset_par(tmp_meas_bt, par_name = "a", dim_meas)
    q_post_me_KK_bt <- get_subset_par(tmp_meas_bt, par_name = "q", dim_meas)
    m_post_me_KK_bt <- get_subset_par(tmp_meas_bt, par_name = "mu", dim_meas)

    mu_info_KK[[regs_to_use[kk]]]   <- compute_mu_info(m_post_me_KK_bt)
    gini_info_KK[[regs_to_use[kk]]] <- compute_gini_info(a_post_me_KK_bt, q_post_me_KK_bt)

    post_marg_kk <- abind::adrop(par_post_estim_fd[, , , kk, drop = FALSE], dim_last_marg)
    tmp_marg_bt  <- f_bt_pars(post_marg_kk, centr_vals, scale_vals)

    a_post_me_fd_bt <- get_subset_par(tmp_marg_bt, par_name = "a", dim_marg)
    q_post_me_fd_bt <- get_subset_par(tmp_marg_bt, par_name = "q", dim_marg)
    m_post_me_fd_bt <- get_subset_par(tmp_marg_bt, par_name = "mu", dim_marg)

    # mu_info_fd_TT[[regs_to_use[kk]]]   <- compute_mu_info(m_post_me_fd_bt) - mu_info_TT_2
    # gini_info_fd_TT[[regs_to_use[kk]]] <- compute_gini_info(a_post_me_fd_bt, q_post_me_fd_bt) - gini_info_TT_2
    # browser()
    mu_info_fd_TT[[regs_to_use[kk]]]   <- compute_mu_info_fd(m_post_me_fd_bt,
                                                             m_post_bt,
                                                             0.8)
    gini_info_fd_TT[[regs_to_use[kk]]] <- compute_gini_info_fd(a_post_me_fd_bt,
                                                               a_post_bt,
                                                               q_post_me_fd_bt,
                                                               q_post_bt,
                                                               0.8)
  }
  return(list(mu_info_TT = mu_info_TT,
              gini_info_TT = gini_info_TT,
              mu_info_fd_TT = mu_info_fd_TT,
              gini_info_fd_TT = gini_info_fd_TT,
              mu_info_KK = mu_info_KK,
              gini_info_KK = gini_info_KK))
}
f_bt_pars <- function(post_estim, values_centering, values_scaling) {
  MM <- dim(post_estim)[3]
  for (mm in seq_len(MM)) {
    post_estim[, , mm] <- post_estim[, , mm] * values_scaling + values_centering
    progress_any(mm, MM)
  }
  return(exp(post_estim))
}
f_bt_marg <- function(post_me, values_centering, values_scaling) {
  MM <- dim(post_me)[3]
  KK <- dim(post_me)[4]
  for (mm in seq_len(MM)) {
    for (kk in seq_len(KK)) {
      post_me[, , mm, kk] <- (post_me[, , mm, kk] * values_scaling +
                                   values_centering)
    }
    progress_any(mm, MM)
  }
  return(exp(post_me))
}
f_bt_meas <- function(post_me, values_centering, values_scaling) {
  GG <- dim(post_me)[3]
  MM <- dim(post_me)[4]
  for (mm in seq_len(MM)) {
    for (gg in seq_len(GG)) {
      post_me[, , gg, mm] <- (post_me[, , gg, mm] * values_scaling +
                                values_centering)
    }
    progress_any(mm, MM)
  }
  return(exp(post_me))
}