# compute_gini_differential_info <- function(out_a, out_q, ki_prob = 0.95) {
#   NN <- dim(out_a)[1]
#   TT <- dim(out_a)[2]
#   if (length(dim(out_a)) > 3) {
#     GG <- dim(out_a)[3]
#   } else {
#     GG <- NULL
#   }
#
#   out_gini <- compute_gini_differential(out_a, out_q, NN, TT, GG, ki_prob = ki_prob)
#   return(out_gini)
# }
# compute_gini_differential <- function(out_a, out_q, NN, TT, GG, ki_prob) {
#   if (is.null(GG)) {
#     out_gini <- array(0, dim = c(NN, TT, 3))
#     for (nn in seq_len(NN)) {
#       for (tt in seq_len(TT)) {
#         vec_tmp <- compute_gini_differential_sm(out_a[nn, tt, ], out_q[nn, tt, ])
#         vec_tmp <- get_vec_tmp(vec_tmp)
#         out_gini[nn, tt, 1]  <- mean(vec_tmp)
#         out_gini[nn, tt, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
#         if (any(is.na(out_gini))) browser()
#       }
#     }
#     new_rn <- get_name_rows_measure_out(rownames(out_a), "_gini")
#     dimnames(out_gini) <- list(new_rn,
#                                colnames(out_a),
#                                c("mean", "ki_upp", "ki_low"))
#   } else {
#     out_gini <- array(0, dim = c(NN, TT, GG, 3))
#     for (nn in seq_len(NN)) {
#       for (tt in seq_len(TT)) {
#         for (gg in seq_len(GG)) {
#           vec_tmp <- compute_gini_differential_sm(out_a[nn, tt, gg, ], out_q[nn, tt, gg, ])
#           vec_tmp <- get_vec_tmp(vec_tmp)
#           out_gini[nn, tt, gg, 1]  <- mean(vec_tmp)
#           out_gini[nn, tt, gg, c(2, 3)] <- quantile(vec_tmp, probs = ki_prob)
#           if (any(is.na(out_gini))) browser()
#         }
#       }
#     }
#     new_rn <- get_name_rows_measure_out(rownames(out_a), "_gini")
#     dimnames(out_gini) <- list(new_rn,
#                                colnames(out_a),
#                                dimnames(out_a)[[3]],
#                                c("mean", "ki_upp", "ki_low"))
#   }
#   return(out_gini)
# }
# compute_gini_differential_sm <- function(a, q, logarithm = FALSE) {
#
# }