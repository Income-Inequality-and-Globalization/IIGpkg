get_cnt_me_B <- function(out_B,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "b", "mu")) {
  nms_rw <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_facs(names_para, num_pars, NN)
  nms_mm <- get_cnt_names_dimMM(dim(out_B)[3])

  dimnames(out_B) <- list(nms_rw, nms_cl, nms_mm)
  return(out_B)
}
get_cnt_me_f <- function(out_f,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "b", "mu")) {
  nms_rw <- get_cnt_names_facs(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_TT(ncol(out_f))
  nms_mm <- get_cnt_names_dimMM(dim(out_f)[3])

  dimnames(out_f) <- list(nms_rw, nms_cl, nms_mm)
  return(out_f)
}
get_cnt_me_D <- function(out_D,
                         NN = 10,
                         num_pars = 3,
                         names_para = c("a", "b", "mu")) {
  nms_rw <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_mm <- get_cnt_names_dimMM(dim(out_D)[3])

  dimnames(out_D) <- list(nms_rw, nms_cl, nms_mm)
  return(out_D)
}
get_cnt_me_wRegs <- function(out_wRegs,
                             NN = 10,
                             num_pars = 3,
                             names_para = c("a", "b", "mu")) {
  nms_rw <- get_cnt_names_NN(names_para, num_pars, NN)
  nms_cl <- get_cnt_names_TT(ncol(out_wRegs))

  rownames(out_wRegs) <- nms_rw
  colnames(out_wRegs) <- nms_cl
  return(out_wRegs)
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
  seq_TT <- seq(from = 1, to = TT, by = 1)
  paste0("tt_", formatC(seq_TT, digits = 1, format = "d", flag = "0"))
}
get_cnt_seq_NN <- function(num_pars, NN) {
  formatC(rep(1:NN, each = num_pars), digits = 1, format = "d", flag = "0")
}
get_cnt_names_dimMM <- function(MM) {
  seq_MM <- seq(from = 1, to = MM, by = 1)
  paste0("MM_", seq_MM)
}