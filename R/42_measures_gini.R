compute_gini_info <- function(out_a, out_q, ki_prob = 0.95) {
  NN <- dim(out_a)[1]
  TT <- dim(out_a)[2]
  if (length(dim(out_a)) > 3) {
    GG <- dim(out_a)[3]
  } else {
    GG <- NULL
  }

  out_gini <- compute_gini(out_a, out_q, NN, TT, GG, ki_prob = ki_prob)
  return(out_gini)
}
compute_gini <- function(out_a, out_q, NN, TT, GG, ki_prob) {
  if (is.null(GG)) {
    out_gini <- array(0, dim = c(NN, TT, 3))
    for (nn in seq_len(NN)) {
      for (tt in seq_len(TT)) {
        vec_tmp <- compute_gini_sm(out_a[nn, tt, ], out_q[nn, tt, ])
        vec_tmp <- get_vec_tmp(vec_tmp, rnd = 1)
        out_gini[nn, tt, 1]  <- mean(vec_tmp)
        out_gini[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
        if (any(is.na(out_gini))) browser()
      }
    }
    new_rn <- get_name_rows_measure_out(rownames(out_a), "_gini")
    dimnames(out_gini) <- list(new_rn,
                               colnames(out_a),
                               c("mean", "ki_upp", "ki_low"))
  } else {
    out_gini <- array(0, dim = c(NN, TT, GG, 3))
    for (nn in seq_len(NN)) {
      for (tt in seq_len(TT)) {
        for (gg in seq_len(GG)) {
          vec_tmp <- compute_gini_sm(out_a[nn, tt, gg, ], out_q[nn, tt, gg, ])
          vec_tmp <- get_vec_tmp(vec_tmp)
          out_gini[nn, tt, gg, 1]  <- mean(vec_tmp)
          out_gini[nn, tt, gg, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
          if (any(is.na(out_gini))) browser()
        }
      }
    }
    new_rn <- get_name_rows_measure_out(rownames(out_a), "_gini")
    dimnames(out_gini) <- list(new_rn,
                               colnames(out_a),
                               dimnames(out_a)[[3]],
                               c("mean", "ki_upp", "ki_low"))
  }
  return(out_gini)
}
compute_gini_sm <- function(a, q, logarithm = FALSE) {
  num_adj <- 0.99999
  if (length(a) != length(q)) browser()
  rhs_top <- lgamma(q) + lgamma(2 * q - 1 / a)
  rhs_low <- lgamma(q - 1 / a) + lgamma(2 * q)

  out_gini_rhs_log <- rhs_top - rhs_low
  out_gini_rhs <- exp(out_gini_rhs_log)

  if (any(out_gini_rhs >= 1)) {
    # out_gini_rhs <- out_gini_rhs[!out_gini_rhs >= 1]
    out_gini_rhs[which(out_gini_rhs >= 1)] <- num_adj
    # stop("Gini values numerically larger 1.")
  }
  if (any(out_gini_rhs <= 0)) {
    # browser()
    # out_gini_rhs <- out_gini_rhs[!out_gini_rhs <= 0]
    out_gini_rhs[which(out_gini_rhs <= 0)] <- 1 - num_adj
    # stop("Gini values numerically smaller 0.")
  }
  if (any(is.na(out_gini_rhs))) {
    # browser()
    out_gini_rhs <- out_gini_rhs[!is.na(out_gini_rhs)]
    # out_gini_rhs[which(out_gini_rhs <= 0)] <- 1 - num_adj
    # stop("Gini values numerically smaller 0.")
  }
  if (isFALSE(logarithm)) out_gini <- 1 - out_gini_rhs
  if (isTRUE(logarithm)) out_gini <- log(1 - out_gini_rhs)
  if (any(is.na(out_gini)) || any(is.nan(out_gini)) || any(out_gini > 1) || any(out_gini < 0)) {
    browser()
  }
  if (length(out_gini) <= 3) browser()
  return(out_gini)
}