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

  B <- get_cnt_me_B(out_gibbs$B)
  f <- get_cnt_me_f(out_gibbs$f)
  D <- get_cnt_me_D(out_gibbs$D)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg)

  # D_mu <- get_subset_par(D, par_name = "mu",
  #                        SUBS_ROWS = TRUE,
  #                        SUBS_COLS = TRUE)
  # B_a <- get_subset_par(B, par_name = "a", SUBS_ROWS = TRUE, SUBS_COLS = TRUE)
  # B_q <- get_subset_par(B, par_name = "q", SUBS_ROWS = TRUE, SUBS_COLS = TRUE)

  NN_num_para <- dim(out_gibbs$B)[1]
  NN <- NN_num_para / const_num_para
  TT <- dim(out_gibbs$f)[2]
  MM <- dim(out_gibbs$B)[3]

  par_post_estim <- get_cnt_par_post_estim(NN = NN, TT = TT, MM = MM)
  for (mm in seq_len(MM)) {
    par_post_estim[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    progress_any(mm, MM)
  }
  par_post_estim_bt <- f_back_transform(par_post_estim, centr_vals, scale_vals)
  a_post_bt  <- get_subset_par(par_post_estim_bt, par_name = "a")
  q_post_bt  <- get_subset_par(par_post_estim_bt, par_name = "q")
  mu_post_bt <- get_subset_par(par_post_estim_bt, par_name = "mu")

  # browser()
  mu_info <- compute_mu_info(mu_post_bt)
  # mu_me_info <- compute_mu_me_info(mu_post_bt, scale_vals, D_mu)
  gini_info <- compute_gini_info(a_post_bt, q_post_bt)
  # browser()
  return(list(mu_info = mu_info,
              gini_info = gini_info))
}
f_back_transform <- function(post_estim, values_centering, values_scaling) {
  # nr <- dim(post_estim)[1]
  # nc <- dim(post_estim)[2]
  MM <- dim(post_estim)[3]
  for (mm in seq_len(MM)) {
    post_estim[, , mm] <- post_estim[, , mm] * values_scaling + values_centering
  }
  return(exp(post_estim))
}