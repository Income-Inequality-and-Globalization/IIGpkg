compute_mu_info <- function(mu_array, ki_prob = 0.95) {
  NN <- dim(mu_array)[1]
  TT <- dim(mu_array)[2]
  if (length(dim(mu_array)) > 3) {
    GG <- dim(mu_array)[3]
  } else {
    GG <- NULL
  }
  out_mu <- compute_mu(mu_array, NN, TT, GG, ki_prob = ki_prob)
  return(out_mu)
}
compute_mu <- function(mu, NN, TT, GG = NULL, ki_prob = ki_prob) {
    if (is.null(GG)) {
      out_mu <- array(0, dim = c(NN, TT, 3))
      for (nn in seq_len(NN)) {
        for (tt in seq_len(TT)) {
          vec_tmp <- get_vec_tmp(mu[nn, tt, ], rnd = 1)
          out_mu[nn, tt, 1]  <- mean(vec_tmp)
          out_mu[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
          if (any(is.na(out_mu))) browser()
        }
      }
      new_rn <- get_name_rows_measure_out(rownames(mu), "_mu")
      dimnames(out_mu) <- list(new_rn,
                               colnames(mu),
                               c("mean", "ki_upp", "ki_low"))
    } else {
      out_mu <- array(0, dim = c(NN, TT, GG, 3))
      for (nn in seq_len(NN)) {
        for (tt in seq_len(TT)) {
          for (gg in seq_len(GG)) {
            vec_tmp <- get_vec_tmp(mu[nn, tt, gg, ])
            out_mu[nn, tt, gg, 1]  <- mean(vec_tmp)
            out_mu[nn, tt, gg, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
          if (any(is.na(out_mu))) browser()
          }
        }
      }
      new_rn <- get_name_rows_measure_out(rownames(mu), "_mu")
      dimnames(out_mu) <- list(new_rn,
                               colnames(mu),
                               dimnames(mu)[[3]],
                               c("mean", "ki_upp", "ki_low"))
    }
  return(out_mu)
}