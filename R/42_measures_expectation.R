compute_mu_info <- function(mu_array, ki_prob = 0.95) {
  NN <- dim(mu_array)[1]
  TT <- dim(mu_array)[2]

  ki_ints <- compute_ki_upper_lower(ki_prob = ki_prob)
  ki_upp <- ki_ints[1]
  ki_low <- ki_ints[2]
  out_mu <- array(0, dim = c(NN, TT, 3))
  for (nn in seq_len(NN)) {
    for (tt in seq_len(TT)) {
      vec_tmp <- mu_array[nn, tt, ]
      out_mu[nn, tt, 1]  <- mean(vec_tmp)
      out_mu[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = c(ki_upp, ki_low))
    }
  }
  dimnames(out_mu) <- list(rownames(mu_array),
                           colnames(mu_array),
                           c("mean", "ki_upp", "ki_low"))
  return(out_mu)
}