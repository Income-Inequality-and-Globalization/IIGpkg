out_to_take <- readRDS(
  paste0("data/output/empirical/Estimates/",
  "cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7/",
  "mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000/",
  "output-Estimates-part2/",
  "cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000.rds"))


# iter_max_mcmc <- dim(out_to_take$Gibbs2_SM_SA$A)[3]
# 
A_means <- apply(out_to_take$A, MARGIN = c(1, 2), mean)
f_means <- apply(out_to_take$f, MARGIN = c(1, 2), mean)
B_means_source <- apply(out_to_take$B, MARGIN = c(1, 2), mean)

B_means <- B_means_source
B_means[, 1] <- B_means[, 1] + (sign(B_means[, 1]))^2
B_means <- B_means[, 2:ncol(B_means)]
diag(B_means) <- diag(B_means) + 1.0

D_means <- apply(out_to_take$D, MARGIN = c(1, 2), mean)
w_regs <- out_to_take$initials$wReg

num_y <- 3
NN <- 10
TT <- 21
NN_TT <- NN * TT

u_simul <- matrix(0, nrow = NN * num_y, ncol = TT)
V_sq_adj <- IIGpkg:::compute_mat_square_root_V(VCOV_array_country_pd, num_y, NN_TT)
set.seed(123)
for (ii in 1:NN) {
  for (tt in 1:TT) {
    Sigma_tmp <- V_sq_adj[, , ((ii - 1) * TT + tt)] %*% A_means %*% t(V_sq_adj[, , ((ii - 1) * TT + tt)])
    if (is.na(Sigma_tmp[1,1])) {
      u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- rep(NA, times = num_y)
    } else {
      u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- MASS::mvrnorm(
      1,
      rep(0, times = num_y),
      V_sq_adj[, , ((ii - 1) * TT + tt)] %*% A_means %*% t(V_sq_adj[, , ((ii - 1) * TT + tt)]))
    }
  }
}

njfacs <- 0
f_means_simul <- matrix(rnorm(NN_TT * num_y + TT * njfacs, mean = 0, sd = 1),
                        nrow = NN * num_y + njfacs, ncol = TT)
f_means_simul <- t(apply(f_means_simul, 1, cumsum))

y_means <- B_means %*% f_means_simul + D_means %*% w_regs

yObs2 <- y_means + u_simul

# plot_title <- paste0("black: data;   ",
#                      "purple: joint fit;   ",
#                      "red: B*f_jnt:   ",
#                      "orange: B*f_idio;   ",
#                      "green: D*w_regs;   ")
# iter_plot <- 1
# plot_name <- paste0("data/output/", "fitted_plot_num_", iter_plot, ".pdf")
# for (t_taken in 1:30) {
#   if (t_taken %% 9 == 1) {
#     pdf(plot_name) 
#     par(mfrow = c(3, 3))
#   }
#   ylim_min <- min(c(yObs2[t_taken, ],
#                     f_all[1, ] * B_means_sim_out[t_taken, 1],
#                     f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1],
#                     (D_means_simul_out %*% w_regs)[t_taken, ],
#                     f_all[1, ] * B_means_sim_out[t_taken, 1] + f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1] + (D_means_simul_out %*% w_regs)[t_taken, ])
#                   , na.rm = TRUE)
#   ylim_max <- max(c(yObs2[t_taken, ],
#                     f_all[1, ] * B_means_sim_out[t_taken, 1],
#                     f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1],
#                     (D_means_simul_out %*% w_regs)[t_taken, ],
#                     f_all[1, ] * B_means_sim_out[t_taken, 1] + f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1] + (D_means_simul_out %*% w_regs)[t_taken, ])
#                   , na.rm = TRUE)
# 
#   plot(yObs2[t_taken, ], type = "l", ylim = c(ylim_min, ylim_max), lwd = 2,
#        ylab =  paste0("No. ", t_taken, " time series"), xlab = "time")
#   
#   # lines((B_means %*% f_means_simul)[t_taken, ], col = "blue")
#   # lines(f_all[1, ], col = "red", lty = "dashed")
#   lines(f_all[1, ] * B_means_sim_out[t_taken, 1], col = "red")
#   # lines(f_all[t_taken + 1, ], col = "orange", lty = "dashed")
#   lines(f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1], col = "orange")
#   # lines(w_regs[t_taken, ], col = "green", lty = "dashed")
#   lines((D_means_simul_out %*% w_regs)[t_taken, ], col = "green")
# 
#   simul_out_fitted <- f_all[1, ] * B_means_sim_out[t_taken, 1] + f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1] + (D_means_simul_out %*% w_regs)[t_taken, ]
#   # simul_out_fitted_no_common <- f_all[1, ] * B_means_sim_out[t_taken, t_taken + 1] + (D_means_simul_out %*% w_regs)[t_taken, ]
#   lines(simul_out_fitted, col = "purple")
#   # lines(simul_out_fitted_no_common, col = "purple", lty = "dashed")
#   if (t_taken %% 9 == 0) {
#     mtext(plot_title, side = 3, line = -2, outer = TRUE)
#     iter_plot <- iter_plot + 1
#     dev.off() 
#     plot_name <- paste0("data/output/", "fitted_plot_num_", iter_plot, ".pdf")
#   }
# }
# mtext(plot_title, side = 3, line = -2, outer = TRUE)
# dev.off() 
# 
# par(mfrow = c(3, 3))
# plot_title <- paste0("black: true simulated factor;   ",
#                      "red: posterior factor;   ")
# iter_plot <- 1
# plot_name <- paste0("data/output/", "factor_plot_num_", iter_plot, ".pdf")
# for (t_taken in 1:31) {
#   if (t_taken %% 9 == 1) {
#     pdf(plot_name) 
#     par(mfrow = c(3, 3))
#   }
#   ylim_min <- min(c(f_means_simul[t_taken, ], f_all[t_taken, ]))
#   ylim_max <- max(c(f_means_simul[t_taken, ], f_all[t_taken, ]))
#   # t_taken <- 17
#   plot(f_means_simul[t_taken, ], type = "l", ylim = c(ylim_min, ylim_max), lwd = 2,
#        ylab =  paste0("No. ", t_taken, " time series"), xlab = "time")
#   lines(f_all[t_taken, ], col = "red")
#   if (t_taken %% 9 == 0) {
#     mtext(plot_title, side = 3, line = -2, outer = TRUE)
#     iter_plot <- iter_plot + 1
#     dev.off() 
#     plot_name <- paste0("data/output/", "factor_plot_num_", iter_plot, ".pdf")
#   }
# }
# mtext(plot_title, side = 3, line = -2, outer = TRUE)
# dev.off() 
# 
# plot_title <- paste0("black: B*f true simulated idio;   ",
#                      "green: B*f posterior factor idio;   ")
# B_times_f_true <- B_means %*% f_means_simul
# B_times_f_post <- B_means_sim_out %*% f_all
# iter_plot <- 1
# plot_name <- paste0("data/output/", "Bf_idio_plot_num_", iter_plot, ".pdf")
# for (t_taken in 1:30) {
#   if (t_taken %% 9 == 1) {
#     pdf(plot_name) 
#     par(mfrow = c(3, 3))
#   }
# 
#     ylim_min <- min(c(f_means_simul[t_taken + 1, ] * B_means[t_taken, t_taken + 1],
#                       f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1]))
#     ylim_max <- max(c(f_means_simul[t_taken + 1, ] * B_means[t_taken, t_taken + 1],
#                       f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1]))
#     plot(f_means_simul[t_taken + 1, ] * B_means[t_taken, t_taken + 1],
#          col = "black", type = "l", ylim = c(ylim_min, ylim_max),
#          ylab = paste0("No. ", t_taken, " time series"))
#     lines(f_all[t_taken + 1, ] * B_means_sim_out[t_taken, t_taken + 1],
#           col = "green")
#   # }
#   if (t_taken %% 9 == 0) {
#     mtext(plot_title, side = 3, line = -2, outer = TRUE)
#     iter_plot <- iter_plot + 1
#     dev.off() 
#     plot_name <- paste0("data/output/", "Bf_idio_plot_num_", iter_plot, ".pdf")
#   }
# }
# mtext(plot_title, side = 3, line = -2, outer = TRUE)
# dev.off() 
# 
# plot_title <- paste0("black: B*f true simulated common;   ",
#                      "red: B*f posterior factor common;   ")
# B_times_f_true <- B_means %*% f_means_simul
# B_times_f_post <- B_means_sim_out %*% f_all
# iter_plot <- 1
# plot_name <- paste0("data/output/", "Bf_common_plot_num_", iter_plot, ".pdf")
# for (t_taken in 1:30) {
#   if (t_taken %% 9 == 1) {
#     pdf(plot_name) 
#     par(mfrow = c(3, 3))
#   }
#     ylim_min <- min(c(f_means_simul[1, ] * B_means[t_taken, 1], f_all[1, ] * B_means_sim_out[t_taken, 1]))
#     ylim_max <- max(c(f_means_simul[1, ] * B_means[t_taken, 1], f_all[1, ] * B_means_sim_out[t_taken, 1]))
#     plot(f_means_simul[1, ] * B_means[t_taken, 1], col = "black", type = "l",
#          ylim = c(ylim_min, ylim_max), ylab =  paste0("No. ", t_taken, " time series"))
#     lines(f_all[1, ] * B_means_sim_out[t_taken, 1], col = "red",
#           ylim = c(ylim_min, ylim_max))
#   if (t_taken %% 9 == 0) {
#     mtext(plot_title, side = 3, line = -2, outer = TRUE)
#     iter_plot <- iter_plot + 1
#     dev.off() 
#     plot_name <- paste0("data/output/", "Bf_common_plot_num_", iter_plot, ".pdf")
#   }
# }
# mtext(plot_title, side = 3, line = -2, outer = TRUE)
# dev.off() 

itermax <- 200000
burnin <- 100000

spec_01 <-  list(
   prior_list = list(
      hyperpriors = list(
        mu_b0 = 0,
        Sigma_b0 = 100,
        alpha_b0 = 10,
        beta_b0 = 0.1
      )
    ),
    OmegaLoad0Scale = 1000,
    OmegaReg0Scale = 1000,
    nu0 = 7,
    Psi0 = diag(3) * 57
)
out <- IIGpkg::Gibbs2_SM_SA_sampler(
      p_joint = njfacs, #p_joint,
      B_par = 0,
      D_par = 0, #grid[j, 6],
      prior_list = spec_01$prior_list,
      nreg = dim(wRegCombsList[[grid[jj, 5]]])[1]/N,
      OmegaLoad0Scale = spec_01$OmegaLoad0Scale, #grid[jj, 1],
      OmegaReg0Scale = spec_01$OmegaReg0Scale, #grid[jj, 2],
      countryA = FALSE,
      A_diag = 1,
      nu0 = spec_01$nu0, # 7, #grid[jj, 3],  #7, #grid[jj, 3],
      Psi0 = spec_01$Psi0, # diag(3) * 57, # grid[jj, 4], # 27,# grid[jj, 4],
      shape0 = shape0,
      rate0 = rate0,
      yObs = yObs2,
      wRegSpec = grid[jj, 5],
      wReg = wRegCombsList[[grid[jj, 5]]],
      Vhat = VCOV_array_country_pd,
      incObsOld = 100000,
      incObsNew = 100000,
      covScale = 1,
      VdiagEst = VdiagEst,
      VhatDiagScale = VhatDiagScale,
      alpha0 = alpha0,
      beta0 = beta0,
      N = N,
      TT = TT,
      storePath = storePath,
      itermax = itermax,
      scaleA = FALSE,
      diagA = FALSE,
      sampleA = TRUE,
      identification = TRUE,
      type = type)

out_gibbs <- out$Gibbs2_SM_SA
# out_gibbs <- readRDS("~/Dropbox/projects/IIG/IIGpkg/data/output/Estimates/cs1_pj0_Reg7_B1_Om1_D0_OmD1_A1_Psi1000_nu1000_IO5000_HYPER/cs1_pj0_B1_Om_D0_OmD1_A1_Psi1000_nu1000_IO5000.rds")
IIGpkg::save_Gibbs_plots(
  Gibbs = out_gibbs,
  burnin = burnin,
  path = plotPath,
  nameMat = nameMat,
  nameReg = nameRegList[[out_gibbs$initials$wRegSpec]],
  njointfac = out_gibbs$initials$njointfac,
  type = type,
  predictionCI = TRUE,
  onlyY = FALSE,
  hierachPrior = TRUE)
