compute_gini_info <- function(out_a, out_q, ki_prob = 0.95) {
  NN <- dim(out_a)[1]
  TT <- dim(out_a)[2]

  ki_ints <- compute_ki_upper_lower(ki_prob = ki_prob)
  ki_upp <- ki_ints[1]
  ki_low <- ki_ints[2]
  out_gini <- array(0, dim = c(NN, TT, 3))
  for (nn in seq_len(NN)) {
    for (tt in seq_len(TT)) {
      vec_tmp <- compute_gini_sm(out_a[nn, tt, ], out_q[nn, tt, ])
      out_gini[nn, tt, 1]  <- mean(vec_tmp)
      out_gini[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = c(ki_upp, ki_low))
    }
  }
  new_rn <- sub("^a_", "", rownames(out_a))
  dimnames(out_gini) <- list(new_rn,
                             colnames(out_a),
                             c("mean", "ki_upp", "ki_low"))
  return(out_gini)
}
compute_gini_sm <- function(a, q, logarithm = FALSE) {
  rhs_top <- lgamma(q) + lgamma(2 * q - 1 / a)
  rhs_low <- lgamma(q - 1 / a) + lgamma(2 * q)

  out_gini_rhs_log <- rhs_top - rhs_low
  out_gini_rhs <- exp(out_gini_rhs_log)
  if (isFALSE(logarithm)) return(1 - out_gini_rhs)
  if (isTRUE(logarithm)) return(log(1 - out_gini_rhs))
}