library(IIGpkg)
load("./data/output/simulations/input_simul_study.RData")
A_means <- input_simul_study$A_means
B_means <- input_simul_study$B_means
D_means <- input_simul_study$D_means
w_regs  <- input_simul_study$w_regs

B_means[, 1] <- B_means[, 1] + (sign(B_means[, 1]))^2
njfacs <- 1
B_means[1, 1] <- 1
diag(B_means[, (1 + njfacs):ncol(B_means)]) <- 1
# B_means <- B_means[, 2:ncol(B_means)]
# diag(B_means) <- diag(B_means) + 1

VCOV_array_country_pd <- input_simul_study$VCOV_array_country_pd
num_y <- 1
NN <- 10
TT <- 21
NN_TT <- NN * TT
V_sq_adj <- IIGpkg:::compute_mat_square_root_V(VCOV_array_country_pd, 3, NN_TT)

NN <- 10
TT <- 100
NN_TT <- TT * NN
NN_num_y <- NN * num_y
nidiofac <- NN * num_y
u_simul    <- matrix(0, nrow = NN * num_y, ncol = TT)
Sigma_true <- array(0, c(num_y, num_y, NN_TT))
V_sq_adj2  <- get_V_sq_adj2(V_sq_adj, NN_TT, num_y, 120)

set.seed(123)
u_simul <- get_u_simul(V_sq_adj2, A_means, num_y, TT, NN)
B_means_used  <- B_means[1:(num_y * NN), 1:(num_y * NN + njfacs)]

f_simul <- generate_f_random_walk(NN_TT, num_y, TT, njfacs)
yObs    <- generate_y_simul(u_simul, B_means_used, f_simul, D = NULL, wRegs = NULL)
plot_data_ml_simul(yObs, cbind(rep(0, times = 10), f_simul[, 2:(TT - 1)]), num_y, NN, njfacs)

Phi   <- diag(NN_num_y + njfacs)
Q     <- diag(NN_num_y + njfacs)
ccc   <- rep(0, NN_num_y)
initX <- rep(0, NN_num_y + njfacs)
initP <- diag(0, NN_num_y + njfacs)
initU <- 0
wReg  <- matrix(rep(1, times = TT), nrow = 1)
VhatArray_A <- array(
      apply(V_sq_adj2[1:num_y, 1:num_y, , drop = FALSE], 3,
            function(X) {
              X %*% as.matrix(A_means[1:num_y, 1:num_y]) %*% t(X)
            }
            ),
      c(num_y, num_y, NN_TT)
    )
VhatArrayBdiagByTime <- IIGpkg::bdiagByTime(VhatArray_A, num_y, NN, TT, NN_num_y)

B_used_tmp <- B_means_used
A_used_tmp <- as.matrix(A_means[1:num_y, 1:num_y])
id_A_mat_optim <- c(diag(matrix(1:num_y^2, nrow = num_y)), which(lower.tri(matrix(1:num_y^2, nrow = num_y))))
num_a_optim <- length(id_A_mat_optim)
Q_used <- Q

true_vals <- c(B_means_used[, 1],
               diag(B_means_used[, (1 + njfacs):ncol(B_means_used)]),
               A_means[id_A_mat_optim])
start_val <- true_vals * 3.5
arg_list_used <- list(yObs = yObs,
                      wReg = wReg,
                      B_used = B_used_tmp,
                      A_used = A_used_tmp,
                      Q_used = Q_used,
                      V_sq_adj = V_sq_adj2,
                      Phi = Phi,
                      ccc = ccc,
                      initX = initX,
                      initU = initU,
                      initP = initP,
                      id_A_mat_optim = id_A_mat_optim,
                      NN = NN,
                      TT = TT,
                      num_y = num_y,
                      NN_TT = NN_TT,
                      njfacs = njfacs)

ml_out_01 <- ml_run_identification_01(true_vals, start_val, HESSIAN = FALSE,
                                      arg_list = arg_list_used)

true_vals_joint_fac <- true_vals[seq_len(NN_num_y)]
optm_vals_joint_fac <- ml_out_01$out_optim$optim_vals[seq_len(NN_num_y)]
par(mfrow = c(3, 3))
gridlength <- 100
grid_width <- 2
for (jj in seq_len(NN_num_y)) {
  if (jj == 1) {
    lower_taken <- 0.0001
    upper_taken <- 2
  } else {
    lower_taken <- NULL
    upper_taken <- NULL
  }
  grid_true_val <- get_seq_grid_ml(true_vals_joint_fac[jj],
                                   lower = lower_taken,
                                   upper = upper_taken,
                                   settings = list(gridlength = gridlength,
                                                   grid_width = grid_width))
  likeli_grid_vals <- vector("numeric", gridlength)
  x_taken <- true_vals
  for (ii in seq_len(gridlength)) {
    x_taken[jj] <- grid_true_val[ii]
    likeli_grid_vals[ii] <- ml_objective_identification_01(
      x_taken,
      yObs = yObs,
      wReg = wReg,
      B_means_used = B_used_tmp,
      A_means_used = A_used_tmp,
      Q_used = Q_used,
      V_sq_adj = V_sq_adj2,
      Phi = Phi,
      ccc = ccc,
      initX = initX,
      initU = initU,
      initP = initP,
      id_A_mat_optim = id_A_mat_optim,
      NN = NN,
      TT = TT,
      num_y = num_y,
      NN_TT = NN_TT,
      num_jnt_fac = njfacs,
      PLOT = TRUE)
  }
  plot(y = likeli_grid_vals, x = grid_true_val, type = "l",
       ylab = "log-likeli",
       xlab = paste0("B_[", jj, ", 1", "] = ",
                     round(B_means_used[jj, 1], digits = 4)),
       main = paste0("ML max.", round(optm_vals_joint_fac[jj], digits = 4)))
  abline(v = true_vals[jj], col = "green")
  abline(v = grid_true_val[which(likeli_grid_vals == max(likeli_grid_vals))],
         col = "red")
    abline(v = optm_vals_joint_fac[jj],
           col = "blue")
  if (jj %% 9 == 0) {
    mtext("KF-loglike on a grid (MAX - RED) ---- true value: green ---- ML. max: blue",
      side = 3, line = -1.5, outer = TRUE)
  }
}

id_B_idio <- seq(from = NN_num_y + 1, to  = 2 * NN_num_y)
true_vals_idio_fac <- true_vals[id_B_idio]
optm_vals_idio_fac <- ml_out_01$out_optim$optim_vals[id_B_idio]
# par(mfrow = c(3, 3))
gridlength <- 100
grid_width <- 0.99
for (jj in seq_len(NN_num_y)) {
  grid_true_val <- get_seq_grid_ml(true_vals_idio_fac[jj],
                                   lower = 0.0001,
                                   upper = 3,
                                   settings = list(gridlength = gridlength,
                                                   grid_width = grid_width))
  likeli_grid_vals <- vector("numeric", gridlength)
  x_taken <- true_vals
  for (ii in seq_len(gridlength)) {
    x_taken[id_B_idio[jj]] <- grid_true_val[ii]
    likeli_grid_vals[ii] <- ml_objective_identification_01(
      x_taken,
      yObs = yObs,
      wReg = wReg,
      B_means_used = B_used_tmp,
      A_means_used = A_used_tmp,
      Q_used = Q_used,
      V_sq_adj = V_sq_adj2,
      Phi = Phi,
      ccc = ccc,
      initX = initX,
      initU = initU,
      initP = initP,
      id_A_mat_optim = id_A_mat_optim,
      NN = NN,
      TT = TT,
      num_y = num_y,
      NN_TT = NN_TT,
      num_jnt_fac = njfacs,
      PLOT = TRUE)
  }
  plot(y = likeli_grid_vals, x = grid_true_val, type = "l",
       ylab = "log-likeli",
       xlab = paste0("B_[", jj, ", ", jj + 1, "] = ",
                     round(B_means_used[jj, njfacs + jj], digits = 4)),
       main = paste0("ML max.", round(optm_vals_idio_fac[jj], digits = 4)))
  abline(v = true_vals_idio_fac[jj], col = "green")
  abline(v = grid_true_val[which(likeli_grid_vals == max(likeli_grid_vals))],
         col = "red")
    abline(v = optm_vals_idio_fac[jj],
           col = "blue")
  if (jj %% 8 == 0) {
    mtext("KF-loglike on a grid (MAX - RED) ---- true value: green ---- ML. max: blue",
      side = 3, line = -1.5, outer = TRUE)
  }
}

id_A <- seq(from = 2 * NN_num_y + 1, to  = 2 * NN_num_y + length(id_A_mat_optim))
true_vals_A <- true_vals[id_A]
optm_vals_A <- ml_out_01$out_optim$optim_vals[id_A]
# par(mfrow = c(1, 1))
gridlength <- 100
grid_width <- 3
for (jj in seq_len(id_A_mat_optim)) {
  grid_true_val <- get_seq_grid_ml(true_vals_A[jj],
                                   upper = 8, lower = 0)
  grid_true_val[which(grid_true_val == 0)] <- grid_true_val[which(grid_true_val == 0) + 1]
  likeli_grid_vals <- vector("numeric", gridlength)
  x_taken <- true_vals
  for (ii in seq_len(gridlength)) {
    x_taken[id_A[jj]] <- grid_true_val[ii]
    likeli_grid_vals[ii] <- ml_objective_identification_01(
      x_taken,
      yObs = yObs,
      wReg = wReg,
      B_means_used = B_used_tmp,
      A_means_used = A_used_tmp,
      Q_used = Q_used,
      V_sq_adj = V_sq_adj2,
      Phi = Phi,
      ccc = ccc,
      initX = initX,
      initU = initU,
      initP = initP,
      id_A_mat_optim = id_A_mat_optim,
      NN = NN,
      TT = TT,
      num_y = num_y,
      NN_TT = NN_TT,
      num_jnt_fac = njfacs,
      PLOT = TRUE)
  }
  plot(y = likeli_grid_vals, x = grid_true_val, type = "l",
       ylab = "log-likeli",
       xlab = paste0("A_[", jj,  ", ", jj, "] = ",
                     round(A_means[jj], digits = 4)),
       main = paste0("ML max.", round(optm_vals_A[jj], digits = 4)))
  abline(v = true_vals_A[jj], col = "green")
  abline(v = grid_true_val[which(likeli_grid_vals == max(likeli_grid_vals))],
         col = "red")
    abline(v = optm_vals_A[jj],
           col = "blue")
}
mtext("KF-loglike on a grid (MAX - RED) ---- true value: green ---- ML. max: blue",
      side = 3, line = -1.5, outer = TRUE)
# 
# gridlength <- 100
# grid_width <- 0.9
# par(mfrow = c(num_y, 2))
# optim_vals_common <- optim_vals[1:(NN_num_y)]
# for (jj in 1:(num_y * NN)) {
#   b_elem <- jj
#   grid_around_belem <- seq(from = B_means_used[b_elem, njfacs + b_elem] - B_means_used[b_elem, njfacs + b_elem] * grid_width,
#                            to = B_means_used[b_elem, njfacs + b_elem] + B_means_used[b_elem, njfacs + b_elem] * grid_width,
#                            length.out = gridlength)
#   iter <- 1
#   likeli_test <- vector("numeric", length(grid_around_belem))
#   B_means_used_tmp <- B_means_used
#   for (gg in grid_around_belem) {
#     B_means_used_tmp[b_elem, njfacs + b_elem] <- gg
#     likeli_test[iter] <- RcppSMCkalman::kfLGSSM(yObs = yObs,
#                                                 uReg = NULL,
#                                                 wReg = wReg,
#                                                 A = Phi,
#                                                 B = NULL,
#                                                 C = B_means_used_tmp,
#                                                 D = ccc,
#                                                 Q = Q, 
#                                                 R = VhatArrayBdiagByTime,
#                                                 initX = initX,
#                                                 initU = initU,
#                                                 initP = initP, 
#                                                 computeMFD = FALSE,
#                                                 computeMSD = FALSE,
#                                                 computePDD = TRUE,
#                                                 computeLLH = TRUE)$kfLLHout
#     iter <- iter + 1
#   }
#   plot(y = likeli_test, x = grid_around_belem, type = "l",
#        ylab = "log-likeli",
#        xlab = paste0("B_idio[", b_elem, ", ", b_elem, "] = ", 
#                      round(B_means_used[b_elem, njfacs + b_elem], digits = 4)),
#        main = paste0("ML max.", ))
#   abline(v = B_means_used[b_elem, njfacs + b_elem], col = "green")
#   abline(v = grid_around_belem[which(likeli_test == max(likeli_test))],
#          col = "red")
#   abline(v = optim_vals_common[jj],
#          col = "blue")
# }
# mtext("KF-loglike on a grid ---- green: true value -- red: KF grid max. -- blue: ML optim.",
#       side = 3, line = -2, outer = TRUE)
# par(mfrow = c(num_y, 2))
# optim_vals_common <- optim_vals[(NN_num_y + 1):(NN_num_y * 2)]
# for (jj in 1:(NN_num_y)) {
#   b_elem <- jj
#   grid_around_belem <- seq(from = B_means_used[b_elem, 1] - B_means_used[b_elem, 1] * grid_width,
#                            to = B_means_used[b_elem, 1] + B_means_used[b_elem, 1] * grid_width,
#                            length.out = gridlength)
#   iter <- 1
#   likeli_test <- vector("numeric", length(grid_around_belem))
#   B_means_used_tmp <- B_means_used
#   for (gg in grid_around_belem) {
#     B_means_used_tmp[b_elem, 1] <- gg
#     likeli_test[iter] <- RcppSMCkalman::kfLGSSM(yObs = yObs,
#                                                 uReg = NULL,
#                                                 wReg = wReg,
#                                                 A = Phi,
#                                                 B = NULL,
#                                                 C = B_means_used_tmp,
#                                                 D = ccc,
#                                                 Q = Q, 
#                                                 R = VhatArrayBdiagByTime,
#                                                 initX = initX,
#                                                 initU = initU,
#                                                 initP = initP, 
#                                                 computeMFD = FALSE,
#                                                 computeMSD = FALSE,
#                                                 computePDD = TRUE,
#                                                 computeLLH = TRUE)$kfLLHout
#     iter <- iter + 1
#   }
#   plot(y = likeli_test, x = grid_around_belem, type = "l",
#        ylab = "log-likeli",
#        xlab = paste0("B_common[", b_elem, ", 1", "] = ", 
#                      round(B_means_used[b_elem, 1], digits = 4)),
#        main = paste0("ML max.", optim_vals_common[jj]))
#   abline(v = B_means_used[b_elem, 1], col = "green")
#   abline(v = grid_around_belem[which(likeli_test == max(likeli_test))],
#          col = "red")
#     abline(v = optim_vals_common[jj],
#          col = "blue")
# }
# mtext("KF-loglike on a grid ---- green: true value ---- red: KF max.",
#       side = 3, line = -2, outer = TRUE)
# par(mfrow = c(num_y, num_y))
# for (jj in 1:num_y) {
#   for (kk in 1:num_y) {
#     a_elem_jj <- jj
#     a_elem_kk <- kk
#     grid_around_aelem <- seq(from = A_means[a_elem_jj, a_elem_kk] - A_means[a_elem_jj, a_elem_kk] * grid_width,
#                              to = A_means[a_elem_jj, a_elem_kk] + A_means[a_elem_jj, a_elem_kk] * grid_width,
#                              length.out = gridlength)
#     iter <- 1
#     likeli_test <- vector("numeric", length(grid_around_aelem))
#     A_means_used_tmp <- A_means
#     for (gg in grid_around_aelem) {
#       A_means_used_tmp[a_elem_jj, a_elem_kk] <- gg
#       A_means_used_tmp[a_elem_kk, a_elem_jj] <- gg
#       VhatArray_A <- array(
#         apply(V_sq_adj2[1:num_y, 1:num_y, , drop = FALSE], 3, 
#               function(X) {X %*% as.matrix(A_means_used_tmp[1:num_y, 1:num_y]) %*% t(X)}),
#         c(num_y, num_y, NN_TT)
#       )
#       VhatArrayBdiagByTime <- IIGpkg::bdiagByTime(VhatArray_A, num_y, NN, TT, NN_num_y)
#       likeli_test[iter] <- RcppSMCkalman::kfLGSSM(yObs = yObs,
#                                                   uReg = NULL,
#                                                   wReg = wReg,
#                                                   A = Phi,
#                                                   B = NULL,
#                                                   C = B_means_used,
#                                                   D = ccc,
#                                                   Q = Q, 
#                                                   R = VhatArrayBdiagByTime,
#                                                   initX = initX,
#                                                   initU = initU,
#                                                   initP = initP, 
#                                                   computeMFD = FALSE,
#                                                   computeMSD = FALSE,
#                                                   computePDD = TRUE,
#                                                   computeLLH = TRUE)$kfLLHout
#       iter <- iter + 1
#     }
#     plot(y = likeli_test, x = grid_around_aelem, type = "l",
#          ylab = "log-likeli",
#          xlab = paste0("A_adjustment[", a_elem_jj, ", ", a_elem_kk, "] = ", 
#                        round(A_means[a_elem_jj, a_elem_kk], digits = 4)))
#     abline(v = A_means[a_elem_jj, a_elem_kk], col = "green")
#     abline(v = grid_around_aelem[which(likeli_test == max(likeli_test))],
#            col = "red")
#   }
# }
# mtext("KF-loglike on a grid ---- green: true value ---- red: KF max.",
#       side = 3, line = -2, outer = TRUE)


# lower_bound <- c(0.0001,
#                  rep(-Inf, times = NN*num_y - 1),
#                  rep(0.0001, times = nidiofac),
#                  rep(0.1, times = num_y),
#                  rep(-Inf, times = num_a_optim - num_y))
# upper_bound <- rep(Inf, times = NN*num_y + nidiofac + num_a_optim)
# check_optim2 <- optimx::optimr(par = start_val,
#                                fn = ml_objective_identification_01,
#                                # gr =  "grnd",
#                                hess = NULL,
#                                method = "L-BFGS-B",
#                                lower = lower_bound,
#                                upper = upper_bound,
#                                hessian = FALSE,
#                                control = list(trace = 3),
#                                yObs = yObs,
#                                wReg = wReg,
#                                B_means_used = B_used_tmp,
#                                A_means_used = A_used_tmp,
#                                Q_used = Q_used,
#                                V_sq_adj = V_sq_adj2,
#                                Phi = Phi,
#                                ccc = ccc,
#                                initX = initX,
#                                initU = initU,
#                                initP = initP,
#                                id_A_mat_optim = id_A_mat_optim,
#                                NN = NN,
#                                TT = TT,
#                                num_y = num_y,
#                                NN_TT = NN_TT,
#                                num_jnt_fac = njfacs)