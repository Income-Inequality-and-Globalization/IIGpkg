#' Computes various a-posterior measures after estimation
#'
#' @param out_gibbs the Gibbs output as an object of class `GibbsOutputIIG`
#' @param transformation_infos the centering and standardizing vectors for the
#'    `y-measurement` data as a vector for all countries (long index) and all
#'    parameters (short index) e.g. a vector of length `30` (`3 x 10`, for ten
#'    countries)
#' @param const_num_para numeric value giving the number of parameters in the
#'    model
#'
#' @return a named list of computed measures
#' @export
generate_measures_me <- function(out_gibbs,
                                 transformation_infos,
                                 const_num_para = 3) {
  if (missing(transformation_infos)) stop("Arg. missing.")
  validate_GibbsOutputIIG(out_gibbs)

  scale_vals  <- transformation_infos$y_settings$scaling
  centr_vals  <- transformation_infos$y_settings$centering
  names_regs  <- transformation_infos$x_settings$names_regs
  cut_off_num <- transformation_infos$x_settings$cut_off_num

  NN_num_par <- dim(out_gibbs$B)[1]
  NN <- NN_num_par / const_num_para
  TT <- dim(out_gibbs$f)[2]
  MM <- dim(out_gibbs$B)[3]
  KK <- length(names_regs)
  GG <- transformation_infos$x_settings$grid_length

  B     <- get_cnt_me_B(out_gibbs$B)
  f     <- get_cnt_me_f(out_gibbs$f)
  D     <- get_cnt_me_D(out_gibbs$D, names_regs)
  D_unc <- get_cnt_me_D_uncertainty(D, names_regs, KK, MM)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg, names_regs)
  WR    <- get_cnt_me_WR(wRegs, names_regs, GG, cut_off_num)
  
  B_m <- apply(B, c(1, 2), mean)
  f_m <- apply(f, c(1, 2), mean)

  par_post_estim <- get_cnt_par_post_estim(NN = NN, TT = TT, MM = MM)
  post_me <- get_cnt_par_me_estim(NN, TT, GG, MM, names_regs)
  for (mm in seq_len(MM)) {
    par_post_estim[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    for (kk in seq_len(KK)) {
      for (gg in seq_len(GG)) {
        post_me[, , gg,  mm, kk] <- B_m %*% f_m + D_unc[ , , mm, kk] %*% WR[, , gg, kk]
      }
    }
    progress_any(mm, MM)
  }
  par_post_estim_bt <- f_bt_par(par_post_estim, centr_vals, scale_vals)

  dim_out <- length(dim(par_post_estim_bt))

  a_post_bt    <- get_subset_par(par_post_estim_bt, par_name = "a", dim_out)
  q_post_bt    <- get_subset_par(par_post_estim_bt, par_name = "q", dim_out)
  mu_post_bt   <- get_subset_par(par_post_estim_bt, par_name = "mu", dim_out)
  mu_info_TT   <- compute_mu_info(mu_post_bt)
  gini_info_TT <- compute_gini_info(a_post_bt, q_post_bt)

  out_mu_info_KK <- vector("list", length = KK)
  names(out_mu_info_KK) <- names_regs
  out_gini_info_KK <- vector("list", length = KK)
  names(out_gini_info_KK) <- names_regs

  for (kk in 1:KK) {
    post_me_kk <- abind::adrop(post_me[, , , , kk, drop = FALSE], 5)
    tmp_me_bt <- f_bt_me(post_me_kk, centr_vals, scale_vals)

    dim_out   <- length(dim(tmp_me_bt))

    a_post_bt <- get_subset_par(tmp_me_bt, par_name = "a", dim_out)
    q_post_bt <- get_subset_par(tmp_me_bt, par_name = "q", dim_out)
    m_post_bt <- get_subset_par(tmp_me_bt, par_name = "mu", dim_out)

    out_mu_info_KK[[kk]]   <- compute_mu_info(m_post_bt)
    out_gini_info_KK[[kk]] <- compute_gini_info(a_post_bt, q_post_bt)
  }
  return(list(mu_info_TT = mu_info_TT,
              gini_info_TT = gini_info_TT,
              mu_info_KK = out_mu_info_KK,
              gini_info_KK = out_gini_info_KK,
              regressor_grid_transformed = WR,
              regressor_grid_original = NULL))
}
f_bt_par <- function(post_estim, values_centering, values_scaling) {
  MM <- dim(post_estim)[3]
  for (mm in seq_len(MM)) {
    post_estim[, , mm] <- post_estim[, , mm] * values_scaling + values_centering
  }
  return(exp(post_estim))
}
f_bt_me <- function(post_me, values_centering, values_scaling) {
  GG <- dim(post_me)[3]
  MM <- dim(post_me)[4]
  for (mm in seq_len(MM)) {
    for (gg in seq_len(GG)) {
      post_me[, , gg, mm] <- post_me[, , gg, mm] * values_scaling + values_centering
    }
  }
  return(exp(post_me))
}