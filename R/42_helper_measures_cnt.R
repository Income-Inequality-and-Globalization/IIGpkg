get_cnt_par_post_estim <- function(NN,
                                   num_pars = 3,
                                   names_para = c("a", "q", "mu"),
                                   TT, MM) {
  out_post_estim <- array(0, dim = c(NN * num_pars, TT, MM))
  dimnames(out_post_estim) <- list(get_cnt_names_NN(names_para, num_pars, NN),
                                   get_cnt_names_TT(TT),
                                   get_cnt_names_dimMM(MM))
  return(out_post_estim)
}
get_cnt_par_me_estim <- function(NN,
                                 TT,
                                 GG,
                                 MM,
                                 names_regs = NULL,
                                 num_pars = 3,
                                 names_para = c("a", "q", "mu")) {
  if (is.null(names_regs)) stop("Arg. 'names_regs' missing.")
  KK <- length(names_regs)
  out_post_estim <- array(0, dim = c(NN * num_pars, TT, GG, MM, KK))
  dimnames(out_post_estim) <- list(get_cnt_names_NN(names_para, num_pars, NN),
                                   get_cnt_names_TT(TT),
                                   get_cnt_names_GG(GG),
                                   get_cnt_names_dimMM(MM),
                                   names_regs)
  return(out_post_estim)
}
get_cnt_me_B <- function(out_B,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "q", "mu")) {
  nms_rw <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_facs(names_para, num_pars, NN)
  nms_mm <- get_cnt_names_dimMM(dim(out_B)[3])

  dimnames(out_B) <- list(nms_rw, nms_cl, nms_mm)
  return(out_B)
}
get_cnt_me_f <- function(out_f,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "q", "mu")) {
  nms_rw <- get_cnt_names_facs(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_TT(ncol(out_f))
  nms_mm <- get_cnt_names_dimMM(dim(out_f)[3])

  dimnames(out_f) <- list(nms_rw, nms_cl, nms_mm)
  return(out_f)
}
get_cnt_me_D <- function(out_D,
                         names_regs,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "q", "mu")) {
  nms_rw <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_regs(names_regs, num_pars, NN)
  nms_mm <- get_cnt_names_dimMM(dim(out_D)[3])

  dimnames(out_D) <- list(nms_rw, nms_cl, nms_mm)
  return(out_D)
}
get_cnt_me_wRegs <- function(out_wRegs,
                             names_regs = NULL,
                             NN = 10,
                             num_pars = 3,
                             names_para = c("a", "q", "mu")) {
  nms_rw <- get_cnt_names_regs(names_regs, num_pars, NN)
  nms_cl <- get_cnt_names_TT(ncol(out_wRegs))
  rownames(out_wRegs) <- nms_rw
  colnames(out_wRegs) <- nms_cl
  return(out_wRegs)
}
get_cnt_me_WR  <- function(wRegs, names_regs = NULL) {
  if (is.null(names_regs)) stop("Arg. 'names_regs' missing.")
  KK <- length(names_regs)
  TT <- ncol(wRegs)
  NN_num_par <- nrow(wRegs)
  # get mean values for each row
  wRegs_means <- apply(wRegs, 1, mean)
  # fill output array with mean values
  out_WR <- array(wRegs_means, dim = c(NN_num_par, TT, KK))
  # fill output array with sorted values at specific rows to get a baseline grid
  for (kk in seq_len(KK)) {
    id_row_kk <- which(grepl(names_regs[kk], rownames(wRegs)))
    for (row_kk in id_row_kk) {
      out_WR[row_kk, , kk] <- sort(wRegs[row_kk, ])
    }
  }
  dimnames(out_WR) <- list(rownames(wRegs), colnames(wRegs), names_regs)
  return(out_WR)
}
get_cnt_names_facs <- function(nms_pars, num_pars, NN) {
  seq_NN <- get_cnt_seq_NN(num_pars, NN)
  c("jF", paste0(paste0("iF", "_", nms_pars), paste0("_NN_", seq_NN)))
}
get_cnt_names_NN <- function(nms_pars, num_pars, NN) {
  seq_NN <- get_cnt_seq_NN(num_pars, NN)
  paste0(nms_pars, "_", "NN_", seq_NN)
}
get_cnt_names_TT <- function(TT) {
  # seq_TT <- seq(from = 1, to = TT, by = 1)
  # paste0("tt_", formatC(seq_TT, digits = 1, format = "d", flag = "0"))
  seq(from = 2000, to = 2000 + TT - 1, by = 1)
}
get_cnt_seq_NN <- function(num_pars, NN) {
  formatC(rep(1:NN, each = num_pars), digits = 1, format = "d", flag = "0")
}
get_cnt_names_dimMM <- function(MM) {
  seq_MM <- seq(from = 1, to = MM, by = 1)
  paste0("MM_", seq_MM)
}
get_cnt_names_regs <- function(names_regs, num_pars, NN) {
  seq_NN <- get_cnt_seq_NN(num_pars, NN)
  paste0("NN_", seq_NN, "_", names_regs)
}