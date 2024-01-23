generate_measures_me <- function(out_gibbs) {
  validate_GibbsOutputIIG(out_gibbs)
  
  B <- get_cnt_me_B(out_gibbs$B)
  f <- get_cnt_me_f(out_gibbs$f)
  D <- get_cnt_me_D(out_gibbs$D)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg)
  
  NN_num_para <- dim(out_test$B)[1]
  TT <- dim(out_test$f)[2]
  MM <- dim(out_test$B)[3]

  par_post_estim <- array(0, dim = c(NN_num_para, TT, MM))
  for (mm in seq_len(MM)) {
    par_post_estim[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    progress_any(mm, MM)
  }
  par_post_estim <- exp(par_post_estim)
  a_post  <- get_par_post_estim(par_post_estim, name_par = "q")
  q_post  <- get_par_post_estim(par_post_estim, name_par = "q")
  mu_post <- get_par_post_estim(par_post_estim, name_par = "mu")
  return(list(B = B, f = f, D = D, wRegs = wRegs,
              par_post_estim = par_post_estim))
}
