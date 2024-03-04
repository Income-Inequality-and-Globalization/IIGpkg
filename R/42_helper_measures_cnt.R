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
  nms_cl <- get_cnt_names_regs_NN(names_regs, num_pars, NN)
  nms_mm <- get_cnt_names_dimMM(dim(out_D)[3])

  dimnames(out_D) <- list(nms_rw, nms_cl, nms_mm)
  return(out_D)
}
get_cnt_me_D_uncertainty <- function(D_mat, names_regs, KK, MM) {
  nrows_tkn  <- nrow(D_mat)
  ncols_tkn  <- ncol(D_mat)
  num_MM_tkn <- dim(D_mat)[[3]]
  D_avg <- array(apply(D_mat, c(1, 2), mean), dim = c(nrows_tkn, ncols_tkn))
  D_out <- array(0, dim = c(dim(D_avg), num_MM_tkn, KK))
  for (kk in seq_len(KK)) {
    id_cols_unc <- get_col_reg_name(colnames(D_mat), names_regs[kk])
    for (mm in seq_len(MM)) {
      D_out[, , mm, kk] <- get_D_uncertain(D_mat[, , mm], D_avg, id_cols_unc)
    }
  }
  dn_from_D_mat <- dimnames(D_mat)
  dimnames(D_out) <- c(dn_from_D_mat, list(get_cnt_names_regs(names_regs, KK)))
  return(D_out)
}
get_cnt_me_wRegs <- function(out_wRegs,
                             names_regs = NULL,
                             NN = 10,
                             num_pars = 3,
                             names_para = c("a", "q", "mu")) {
  nms_rw <- get_cnt_names_regs_NN(names_regs, num_pars, NN)
  nms_cl <- get_cnt_names_TT(ncol(out_wRegs))
  rownames(out_wRegs) <- nms_rw
  colnames(out_wRegs) <- nms_cl
  return(out_wRegs)
}
get_cnt_me_WR  <- function(wRegs, names_regs = NULL, grid_length, cutoff_num) {
  if (is.null(names_regs)) stop("Arg. 'names_regs' missing.")
  KK <- length(names_regs)
  GG <- grid_length
  NN_num_par <- nrow(wRegs)
  TT <- ncol(wRegs)
  # out_WR <- array(wRegs, dim = c(NN_num_par, TT, GG, KK))
  # get mean values for each row
  wRegs_means <- apply(wRegs, 1, mean)
  # fill output array with mean values
  out_WR <- array(wRegs_means, dim = c(NN_num_par, TT, GG, KK))
  # fill output array with sorted values at specific rows to get a baseline grid
  for (kk in seq_len(KK)) {
    id_row_kk <- get_col_reg_name(rownames(wRegs), names_regs[kk])
    for (row_kk in id_row_kk) {
      out_WR[row_kk, , , kk] <- get_single_reg_grid(
        wRegs[row_kk, ],
        grid_length,
        cutoff_num)
    }
  }
  dimnames(out_WR) <- list(rownames(wRegs),
                           get_cnt_names_TT(TT),
                           get_cnt_names_GG(GG),
                           names_regs)
  return(out_WR)
}
get_single_reg_grid <- function(vals, grid_length, cutoff_num) {
  TT <- length(vals)
  out_reg_grid <- matrix(0, nrow = TT, ncol = grid_length)
  for (tt in seq_len(TT)) {
    out_reg_grid[tt, ] <- get_grid_vals(vals, grid_length, cutoff_num, tt)
  }
  return(out_reg_grid)
}
get_grid_vals <- function(values, grid_length, cutoff, tt) {
  ##############################################################################
  # THIS IS JUST TESTING
  # if (values[tt] == 0) values[tt] <- 0.0001
  # out_grid_vals <- rep(values[tt], times = grid_length)
  # THIS IS JUST TESTING
  ##############################################################################

  # out_grid_vals <- get_grid_vals_small_large_cutoff(values,
  #                                                   grid_length,
  #                                                   cutoff)
  # out_grid_vals <- get_grid_vals_small_large_nonzero(values,
  #                                                    grid_length,
  #                                                    cutoff)
  out_grid_vals <- get_grid_vals_remove_intervall(values, grid_length, cutoff)
  # out_grid_vals <- get_grid_vals_around_mean(values, grid_length, tt)
  return(out_grid_vals)
}
get_grid_vals_around_mean <- function(values, grid_length, tt) {
  cntr_val <- values[tt]
  around_upp_low <- c(cntr_val + c(-0.1, 0.1) * sd(values))
  out_grid_vals <- seq(from = around_upp_low[1],
                       to = around_upp_low[2],
                       length.out = grid_length)
  if (any(out_grid_vals == 0)) stop("Grid values should not be zero.")
  return(out_grid_vals)
}
get_grid_vals_remove_intervall <- function(x, grid_length, cutoff) {
  grid_vals <- sort(x)
  grid_vals <- get_vals_nozero_intervall(grid_vals)
  # grid_vals <- get_vals_cutoff(grid_vals, cutoff)
  if (any(grid_vals == 0)) stop("Grid values should not be zero.")
  grid_vals_neg <- grid_vals[grid_vals < 0]
  grid_vals_pos <- grid_vals[grid_vals > 0]
  num_neg <- length(grid_vals_neg)
  num_pos <- length(grid_vals_pos)
  if (num_neg >= 2) {
    len_out   <- floor(num_neg / length(grid_vals) * grid_length)
    min_start <- min(grid_vals_neg)
    max_end   <- max(grid_vals_neg)
    if (min_start == max_end) max_end <- max_end - 0.1
    grid_vals_neg <- seq(from = min_start, to = max_end, length.out = len_out)
    if (length(grid_vals_neg) != len_out) {
      browser()
    }
  } else {
    len_out <- num_neg
  }
  if (num_pos >= 1) {
    len_out <- grid_length - len_out
    # if (num_pos == 1 && len_out > 2) browser()
    min_start <- min(grid_vals_pos)
    max_end   <- max(grid_vals_pos)
    if (min_start == max_end) max_end <- max_end + 0.1 * (len_out - 1)
    grid_vals_pos <- seq(from = min_start, to = max_end, length.out = len_out)
  } else {
    len_out <- grid_length - len_out
    if (len_out != num_pos) {
      grid_vals_pos <- c(grid_vals_pos, grid_vals_pos + 0.1)
    }
  }
  out_grid_vals <- c(grid_vals_neg, grid_vals_pos)
  if (length(out_grid_vals) != grid_length) browser()
  return(unname(out_grid_vals))
}
get_vals_nozero_intervall <- function(vals_to_trim, trim_config = 2) {
  id_to_trim <- which(vals_to_trim == 0)
  if (length(id_to_trim) == 0) {
    id_to_trim <- which(min(vals_to_trim))
  }
  if (length(id_to_trim) > 1) id_to_trim <- id_to_trim[1]

  b_low <- vals_to_trim[id_to_trim] - sd(vals_to_trim) / trim_config
  b_upp <- vals_to_trim[id_to_trim] + sd(vals_to_trim) / trim_config
  out_vals <- vals_to_trim[!(vals_to_trim >= b_low & vals_to_trim <= b_upp)]

  return(out_vals)
}
get_grid_vals_small_large_nonzero <- function(x, grid_length, cutoff) {
  stopifnot(`Something is wrong; there must be zero regs vals.` = any(x == 0))
  out_grid_vals <- sort(x[x != 0])
  if (cutoff != 0 || is.null(cutoff)) {
    out_grid_vals <- out_grid_vals[-c(1:cutoff)]
  }
  if (grid_length != length(out_grid_vals)) {
    out_grid_vals <- append_to_vector(out_grid_vals,
                                      mean(diff(out_grid_vals)),
                                      grid_length,
                                      from = "top")
  }
  if (any(out_grid_vals == 0)) stop("Grid values should not be zero.")
  return(out_grid_vals)
}
append_to_vector <- function(vec, increment, final_length, from = "bottom") {
  if (length(vec) >= final_length) {
    warning("Initial vector is already of the required length or longer.")
    return(vec)
  }
  while (length(vec) < final_length) {
    if (length(vec) %% 2 == 0 && from == "both") {
      # Add to the bottom
      vec <- c(vec[1] - increment, vec)
    } else if (from == "both") {
      # Add to the top
      vec <- c(vec, tail(vec, 1) + increment)
    } else if (from == "bottom") {
      # Add to the bottom
      vec <- c(vec[1] - increment, vec)
    } else if (from == "top") {
      # Add to the top
      vec <- c(vec, tail(vec, 1) + increment)
    }
  }
  return(unname(vec))
}
get_grid_vals_small_large_cutoff <- function(x, grid_length, cutoff) {
  x <- sort(x)
  x <- get_vals_cutoff(x, cutoff)
  v_max <- max(x)
  v_min <- min(x)
  out_grid_vals <- seq(from = v_min, to = v_max, length.out = grid_length)
  out_grid_vals <- set_zero_to_closest_nonzero(out_grid_vals)
  if (any(out_grid_vals == 0)) stop("Grid values should not be zero.")
  return(out_grid_vals)
}
get_vals_cutoff <- function(vals_to_cut, cutoff) {
  browser()
  if (cutoff == 0 || is.null(cutoff)) {
    from_low <- 0
    from_upp <- 0
  } else {
    from_low <- seq_len(cutoff)
    from_upp <- length(vals_to_cut):(length(vals_to_cut) - cutoff + 1)
  }
  id_to_use <- setdiff(seq_along(vals_to_cut), c(from_low, from_upp))
  vals_to_use <- vals_to_cut[id_to_use]
  return(vals_to_use)
}
set_zero_to_closest_nonzero <- function(vec) {
  id_replace <- vec == 0
  vec_check <- vec[vec != 0]
  id_cl_friend <- which(abs(0 - vec_check) == min(abs(0 - vec_check)))
  vec[id_replace] <- vec_check[id_cl_friend] / 2
  return(vec)
}
get_cnt_info_KK <- function(KK, names_regs) {
  out_mu_info_KK <- vector("list", length = KK)
  names(out_mu_info_KK) <- names_regs
  return(out_mu_info_KK)
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
  seq(from = 2000, to = 2000 + TT - 1, by = 1)
}
get_cnt_names_GG <- function(GG) {
  seq_GG <- seq(from = 1, to = GG, by = 1)
  paste0("gg_", formatC(seq_GG, digits = 1, format = "d", flag = "0"))
}
get_cnt_seq_NN <- function(num_pars, NN) {
  formatC(rep(1:NN, each = num_pars), digits = 1, format = "d", flag = "0")
}
get_cnt_names_dimMM <- function(MM) {
  seq_MM <- seq(from = 1, to = MM, by = 1)
  paste0("MM_", seq_MM)
}
get_cnt_names_regs_NN <- function(names_regs, num_pars, NN) {
  seq_NN <- get_cnt_seq_NN(num_pars, NN)
  paste0("NN_", seq_NN, "_", names_regs)
}
get_cnt_names_regs <- function(names_regs, KK) {
  paste0("KK_",
         formatC(seq_len(KK), digits = 1, format = "d", flag = "0"),
         "_",
         names_regs)
}
get_col_reg_name <- function(all_names_Dcol, name_reg) {
  which(grepl(name_reg, all_names_Dcol))
}
get_D_uncertain <- function(D, D_avg, id_cu) {
  D_avg[, id_cu] <- D[, id_cu]
  return(D_avg)
}