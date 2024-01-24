generate_measures_me <- function(out_gibbs,
                                 transformation_infos,
                                 const_num_para = 3) {
  if (missing(transformation_infos)) stop("Arg. missing.")
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
  par_post_estim <- f_back_transform(par_post_estim,
                                     transformation_infos$centering,
                                     transformation_infos$scaling)
  a_post  <- get_par_post_estim(par_post_estim, par_name = "a")
  q_post  <- get_par_post_estim(par_post_estim, par_name = "q")
  mu_post <- get_par_post_estim(par_post_estim, par_name = "mu")

  mu_info <- compute_mu_info(mu_post)
  gini_info <- compute_gini_info(a_post, q_post)
  return(list(B = B, f = f, D = D, wRegs = wRegs,
              par_post_estim = par_post_estim,
              mu_info = mu_info, gini_info = gini_info))
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