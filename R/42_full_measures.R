generate_measures_me <- function(out_gibbs, const_num_para = 3) {
  validate_GibbsOutputIIG(out_gibbs)

  B <- get_cnt_me_B(out_gibbs$B)
  f <- get_cnt_me_f(out_gibbs$f)
  D <- get_cnt_me_D(out_gibbs$D)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg)

  NN_num_para <- dim(out_gibbs$B)[1]
  NN <- NN_num_para / const_num_para
  TT <- dim(out_gibbs$f)[2]
  MM <- dim(out_gibbs$B)[3]

  par_post_estim <- get_cnt_par_post_estim(NN = NN, TT = TT, MM = MM)
  for (mm in seq_len(MM)) {
    par_post_estim[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    progress_any(mm, MM)
  }
  par_post_estim <- exp(par_post_estim)
  a_post  <- get_par_post_estim(par_post_estim, par_name = "a")
  q_post  <- get_par_post_estim(par_post_estim, par_name = "q")
  mu_post <- get_par_post_estim(par_post_estim, par_name = "mu")

  mu_info <- compute_mu_info(mu_post)
  browser()
  gini_info <- compute_gini_info(a_post, q_post)
  return(list(B = B, f = f, D = D, wRegs = wRegs,
              par_post_estim = par_post_estim,
              mu_info = mu_info, gini_info = gini_info))
}
