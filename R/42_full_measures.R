generate_measures_me <- function(out_gibbs) {
  validate_GibbsOutputIIG(out_gibbs)
  
  B <- get_cnt_me_B(out_gibbs$B)
  f <- get_cnt_me_f(out_gibbs$f)
  D <- get_cnt_me_D(out_gibbs$D)
  wRegs <- get_cnt_me_wRegs(out_gibbs$initials$wReg)
  
  NN_num_para <- dim(out_test$B)[1]
  TT <- dim(out_test$f)[2]
  MM <- dim(out_test$B)[3]

  out_mu <- array(0, dim = c(NN_num_para, TT, MM))
  for (mm in seq_len(MM)) {
    out_mu[, , mm] <- B[, , mm] %*% f[, , mm] + D[, , mm] %*% wRegs
    progress_any(mm, MM)
  }
  out_mu <- exp(out_mu)
  return(list(B = B, f = f, D = D, wRegs = wRegs, mu = out_mu))
}
