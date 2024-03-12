get_subset_par <- function(out_post_all,
                           par_name,
                           dim_out,
                           SUBS_ROWS = TRUE,
                           SUBS_COLS = FALSE
                           # id_keep_cols = NULL,
                           # id_keep_rows = NULL
                           ) {
  stopifnot(`Invalid argument 'par_name'.` = par_name %in% c("a", "q", "mu"))
  if (SUBS_ROWS) {
    row_names_par <- rownames(out_post_all)
    id_taken <- find_id_par_name_all(par_name, row_names_par)
    # if (!is.null(id_keep_rows)) id_taken <- sort(c(id_taken, id_keep_rows))
    if (dim_out == 3) out_post_all <- out_post_all[id_taken, , , drop = FALSE]
    if (dim_out == 4) out_post_all <- out_post_all[id_taken, , , , drop = FALSE]
  }
  if (SUBS_COLS) {
    col_names_par <- colnames(out_post_all)
    id_taken <- find_id_par_name_all(par_name, col_names_par)
    # if (!is.null(id_keep_cols)) id_taken <- sort(c(id_taken, id_keep_cols))
    if (dim_out == 3) out_post_all <- out_post_all[, id_taken, , drop = FALSE]
    if (dim_out == 4)  out_post_all <- out_post_all[, id_taken, , , drop = FALSE]
  }
  return(out_post_all)
}
find_id_par_name_all <- function(par_name, names_to_search) {
  tmp_rgx <- get_regex_par(par_name)
  which(grepl(tmp_rgx, names_to_search))
}
get_regex_par <- function(par_name) {
  if (par_name == "a") return("a")
  if (par_name == "q") return("q")
  if (par_name == "mu") return("mu")
}
get_names_measures <- function(nm_measure = NULL) {
  nm_measure <- paste0(nm_measure, "_")
  list(mean = paste0(nm_measure, "mean"),
       ki_low = paste0(nm_measure, "ki_low"),
       ki_upp = paste0(nm_measure, "ki_upp"))
}
get_vec_tmp <- function(x, rnd = 1.0) {
  x2 <- sort(x)
  ct  <- floor((length(x) - length(x) * rnd) / 2)
  if (ct != 0) {
    out <- x2[-c(1:ct, (length(x) - ct + 1):length(x))]
  } else {
    out <- x2
  }
  # out <- setdiff(x, c(head(x, n = ct), tail(x, n = ct)))
  if (length(out) == 0) browser()
  return(out)
  return(x)
}
get_name_rows_measure_out <- function(rownames, suffix) {
  clean_strng <- paste0(sub("^a_", "", rownames))
  clean_strng <- paste0(sub("^b_", "", clean_strng))
  clean_strng <- paste0(sub("^q_", "", clean_strng))
  clean_strng <- paste0(sub("^mu_", "", clean_strng))
  paste0(clean_strng, suffix)
}
get_regs_grid_TT <- function(regs_grid_TT_used, names_regs) {
  num_regs <- length(names_regs)
  NN <- dim(regs_grid_TT_used)[1] / num_regs
  TT <- dim(regs_grid_TT_used)[2]
  out_regs_array <- array(0, dim = c(NN, TT, num_regs))
  dimnames(out_regs_array) <- list(
    paste0("NN_", formatC(1:NN, digits = 1, format = "d", flag = "0")),
    colnames(regs_grid_TT_used),
    names_regs)
  for (kk in seq_along(names_regs)) {
    id_regs <- seq(from = kk, to = num_regs * (NN - 1) + kk, by = 3)
    out_regs_array[, , kk] <- regs_grid_TT_used[id_regs, ]
  }
  return(out_regs_array)
}
adjust_pnt <- function(pnt_tkn, plot_data, thr = 0.3) {
  if (is.na(pnt_tkn[["x"]])) return(pnt_tkn)
  abs_diff <- abs(plot_data$x - pnt_tkn[["x"]])
  min_grid_pnt <- min(abs_diff)
  id_closest_pnt <- which(min_grid_pnt == abs_diff)

  CHECK_ADJUST <- 100 * abs(pnt_tkn[["y"]] - plot_data$y[id_closest_pnt]) / plot_data$y[id_closest_pnt] > thr
  if (!is.logical(CHECK_ADJUST)) {browser()}
  if (length(CHECK_ADJUST) == 0) {browser()}

  if (CHECK_ADJUST) {
  pnt_tkn[["y"]] <-  plot_data$y[id_closest_pnt]
  }
  return(pnt_tkn)
}
