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

  scale_vals <- transformation_infos$scaling
  centr_vals <- transformation_infos$centering
  names_regs <- transformation_infos$names_regs
  KK         <- length(names_regs)

  B     <- get_cnt_me_B(out_gibbs$B)
  f     <- get_cnt_me_f(out_gibbs$f)
  D     <- get_cnt_me_D(out_gibbs$D)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg, names_regs)
  WR    <- get_cnt_me_WR(wRegs, names_regs)

  NN_num_par <- dim(out_gibbs$B)[1]
  NN <- NN_num_par / const_num_para
  TT <- dim(out_gibbs$f)[2]
  MM <- dim(out_gibbs$B)[3]

  par_post_estim <- get_cnt_par_post_estim(NN = NN, TT = TT, MM = MM)
  post_me <- get_cnt_par_me_estim(NN, names_regs, TT, MM)
  for (mm in seq_len(MM)) {
    par_post_estim[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    for (kk in seq_len(KK)) {
      post_me[, , mm, kk] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% WR[, , kk]
    }
    progress_any(mm, MM)
  }
  par_post_estim_bt <- f_back_transform(par_post_estim, centr_vals, scale_vals)
  a_post_bt    <- get_subset_par(par_post_estim_bt, par_name = "a")
  q_post_bt    <- get_subset_par(par_post_estim_bt, par_name = "q")
  mu_post_bt   <- get_subset_par(par_post_estim_bt, par_name = "mu")
  mu_info_TT   <- compute_mu_info(mu_post_bt)
  gini_info_TT <- compute_gini_info(a_post_bt, q_post_bt)

  out_mu_info_KK <- vector("list", length = KK)
  names(out_mu_info_KK) <- names_regs
  out_gini_info_KK <- vector("list", length = KK)
  names(out_gini_info_KK) <- names_regs
  for (kk in 1:3) {
    tmp_me_bt <- f_back_transform(post_me[, , , kk], centr_vals, scale_vals)
    a_post_bt <- get_subset_par(tmp_me_bt, par_name = "a")
    q_post_bt <- get_subset_par(tmp_me_bt, par_name = "q")
    m_post_bt <- get_subset_par(tmp_me_bt, par_name = "mu")
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
f_back_transform <- function(post_estim, values_centering, values_scaling) {
  MM <- dim(post_estim)[3]
  for (mm in seq_len(MM)) {
    post_estim[, , mm] <- post_estim[, , mm] * values_scaling + values_centering
  }
  return(exp(post_estim))
}