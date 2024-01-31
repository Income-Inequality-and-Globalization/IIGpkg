compute_mu_info <- function(mu_array, ki_prob = 0.95) {
  NN <- dim(mu_array)[1]
  TT <- dim(mu_array)[2]
  if (length(dim(mu_array)) > 3) {
    GG <- dim(mu_array)[3]
  } else {
    GG <- NULL
  }
  ki_int <- compute_ki_upper_lower(ki_prob = ki_prob)
  ki_upp <- ki_int[1]
  ki_low <- ki_int[2]
  out_mu <- compute_mu(mu_array, NN, TT, GG, ki_upp, ki_low)
  return(out_mu)
}
compute_mu <- function(mu, NN, TT, GG = NULL, ki_upp, ki_low) {
    if (is.null(GG)) {
      out_mu <- array(0, dim = c(NN, TT, 3))
      for (nn in seq_len(NN)) {
        for (tt in seq_len(TT)) {
          vec_tmp <- mu[nn, tt, ]
          out_mu[nn, tt, 1]  <- mean(vec_tmp)
          out_mu[nn, tt, c(2, 3)] <- quantile(vec_tmp,
                                              probs = c(ki_upp, ki_low))
        }
      }
      dimnames(out_mu) <- list(rownames(mu),
                               colnames(mu),
                               c("mean", "ki_upp", "ki_low"))
    } else {
      out_mu <- array(0, dim = c(NN, TT, GG, 3))
      for (gg in seq_len(GG)) {
        for (nn in seq_len(NN)) {
          for (tt in seq_len(TT)) {
            vec_tmp <- mu[nn, tt, gg, ]
            out_mu[nn, tt, gg, 1]  <- mean(vec_tmp)
            out_mu[nn, tt, gg, c(2, 3)] <- quantile(vec_tmp,
                                                    probs = c(ki_upp, ki_low))
          }
        }
      }
      dimnames(out_mu) <- list(rownames(mu),
                               colnames(mu),
                               dimnames(mu)[[3]],
                               c("mean", "ki_upp", "ki_low"))
    }
  return(out_mu)
}
# compute_mu_me_info <- function(mu_array_bt,
#                                scale_values,
#                                loadings,
#                                ki_prob = 0.95) {
#   NN <- dim(mu_array_bt)[1]
#   TT <- dim(mu_array_bt)[2]
# 
#   ki_int <- compute_ki_upper_lower(ki_prob = ki_prob)
#   ki_upp <- ki_int[1]
#   ki_low <- ki_int[2]
#   out_mu <- array(0, dim = c(NN, TT, 3))
#   for (nn in seq_len(NN)) {
#     for (tt in seq_len(TT)) {
#       vec_tmp <- mu_array_bt[nn, tt, ]
#       out_mu[nn, tt, 1]  <- mean(vec_tmp)
#       out_mu[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = c(ki_upp, ki_low))
#     }
#   }
#   dimnames(out_mu) <- list(rownames(mu_array_bt),
#                            colnames(mu_array_bt),
#                            c("mean", "ki_upp", "ki_low"))
#   return(out_mu)
# }
# compute_mu_me <- function(mu_vec, loading_vec, scale_vls) {
#   mu_vec * loading_vec * scale_vls
# }