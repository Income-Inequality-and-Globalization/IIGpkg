library(IIGpkg)
load("./data/output/simulations/input_simul_study.RData")
A_means <- input_simul_study$A_means
B_means <- input_simul_study$B_means
D_means <- input_simul_study$D_means
w_regs  <- input_simul_study$w_regs

B_means[, 1] <- B_means[, 1] + (sign(B_means[, 1]))^2
njfacs <- 1
B_means[1, 1] <- 1
# diag(B_means[, (1 + njfacs):ncol(B_means)]) <- 1
VCOV_array_country_pd <- input_simul_study$VCOV_array_country_pd

num_y <- 3
NN <- 10
TT <- 50
# TT <- 21
NN_TT <- NN * TT

VCOV_array_country_pd <- adjust_vcov_array_country_pd_to_TT(
  VCOV_array_country_pd, 
  TT,
  NN)
w_regs <- adjust_w_regs_to_TT(
  w_regs, 
  TT)


u_simul <- matrix(0, nrow = NN * num_y, ncol = TT)
V_sq_adj <- IIGpkg:::compute_mat_square_root_V(VCOV_array_country_pd, num_y, NN_TT)
Sigma_true <- array(0, c(3, 3, NN_TT))
set.seed(123)
it <- 1
for (ii in 1:NN) {
  for (tt in 1:TT) {
    Sigma_tmp <- V_sq_adj[, , ((ii - 1) * TT + tt)] %*% A_means %*% t(V_sq_adj[, , ((ii - 1) * TT + tt)])
    Sigma_true[,, it] <- Sigma_tmp
    if (is.na(Sigma_tmp[1,1])) {
      u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- rep(NA, times = num_y)
    } else {
      u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- MASS::mvrnorm(
      1,
      rep(0, times = num_y),
      V_sq_adj[, , ((ii - 1) * TT + tt)] %*% A_means %*% t(V_sq_adj[, , ((ii - 1) * TT + tt)]))
    }
    it <- it + 1
  }
}

f_means_simul <- matrix(rnorm(NN_TT * num_y + TT * njfacs, mean = 0, sd = 1),
                        nrow = NN * num_y + njfacs, ncol = TT)
f_means_simul <- t(apply(f_means_simul, 1, cumsum))
y_means <- B_means %*% f_means_simul + D_means %*% w_regs
yObs <- y_means + u_simul

pth_true_vals <- "./kf_test2.RData"
true_vals <- set_true_vals(pth_true_vals,
                           VCOV_array_country_pd,
                           A_means,
                           B_means,
                           D_means,
                           yObs,
                           w_regs,
                           f_means_simul)

# itermax <- 50000
# burnin  <- 10000

itermax <- 10000
burnin  <- 5000

spec_01 <- list(
   prior_list = NULL,
    #  list(
    #   hyperpriors = list(
    #     mu_b0 = 0,
    #     Sigma_b0 = 100,
    #     alpha_b0 = 10,
    #     beta_b0 = 0.1
    #   )
    # ),
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
      nreg = input_simul_study$nreg, # dim(wRegCombsList[[grid[jj, 5]]])[1]/N,
      OmegaLoad0Scale = spec_01$OmegaLoad0Scale, #grid[jj, 1],
      OmegaReg0Scale = spec_01$OmegaReg0Scale, #grid[jj, 2],
      countryA = FALSE,
      A_diag = 1,
      nu0 = spec_01$nu0, # 7, #grid[jj, 3],  #7, #grid[jj, 3],
      Psi0 = spec_01$Psi0, # diag(3) * 57, # grid[jj, 4], # 27,# grid[jj, 4],
      shape0 = input_simul_study$shape0,
      rate0 = input_simul_study$rate0,
      yObs = yObs,
      wRegSpec = 7, #grid[jj, 5],
      wReg = w_regs, # wRegCombsList[[grid[jj, 5]]],
      Vhat = VCOV_array_country_pd, # VCOV_array_country_pd
      incObsOld = 100000,
      incObsNew = 100000,
      covScale = 1,
      VdiagEst = input_simul_study$VdiagEst,
      VhatDiagScale = input_simul_study$VhatDiagScale,
      alpha0 = input_simul_study$alpha0,
      beta0 = input_simul_study$beta0,
      N = NN,
      TT = TT,
      storePath = input_simul_study$storePath,
      itermax = itermax,
      scaleA = FALSE,
      diagA = FALSE,
      sampleA = TRUE,
      identification = TRUE,
      type = input_simul_study$type)

# out$Gibbs2_SM_SA$f <- array(replicate(4, out$Gibbs2_SM_SA$f), dim = c(31, 50, 200000))
# out$Gibbs2_SM_SA$B <- array(replicate(4, out$Gibbs2_SM_SA$B), dim = c(30, 31, 200000))
# out$Gibbs2_SM_SA$D <- array(replicate(4, out$Gibbs2_SM_SA$D), dim = c(30, 30, 200000))
# out$Gibbs2_SM_SA$A <- array(replicate(4, out$Gibbs2_SM_SA$A), dim = c(3, 3, 200000))
# out$Gibbs2_SM_SA$V <- array(replicate(4, out$Gibbs2_SM_SA$V), dim = c(30, 200000))
# out$Gibbs2_SM_SA$BD0STORE <- array(replicate(4, out$Gibbs2_SM_SA$BD0STORE), dim = c(15, 200000))
# out$Gibbs2_SM_SA$Omega0STORE <- array(replicate(4, out$Gibbs2_SM_SA$Omega0STORE), dim = c(15, 200000))
# out$Gibbs2_SM_SA$initials$itermax <- out$Gibbs2_SM_SA$initials$itermax*4
# saveRDS(out$Gibbs2_SM_SA,
#         file = "~/Dropbox/projects/IIG/IIGpkg/data/output/simulations/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7_IO1e+05_TT50_AT_0.1_RESTR-N1/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7_IO1e+05.rds")
# out_gibbs$f <- array(replicate(5, out_gibbs$f[, , 10001:50000]), dim = c(31, 50, 200000))
# out_gibbs$B <- array(replicate(5, out_gibbs$B[, , 10001:50000]), dim = c(30, 31, 200000))
# out_gibbs$D <- array(replicate(5, out_gibbs$D[, , 10001:50000]), dim = c(30, 30, 200000))
# out_gibbs$A <- array(replicate(5, out_gibbs$A[, , 10001:50000]), dim = c(3, 3, 200000))
# out_gibbs$V <- array(replicate(4, out_gibbs$V), dim = c(30, 200000))
# out_gibbs$BD0STORE <- array(replicate(4, out_gibbs$BD0STORE), dim = c(15, 200000))
# out_gibbs$Omega0STORE <- array(replicate(4, out_gibbs$Omega0STORE), dim = c(15, 200000))
out_gibbs <- out$Gibbs2_SM_SA
# out_gibbs <- readRDS(
#   file.path(
#     "/home/iz/Dropbox/projects/IIG/IIGpkg/data",
#     "output/simulations/Estimates",
#     "cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7_IO1e+05_TT50_AT_ZERO_RESTR-N1",
#     "cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7_IO1e+05.rds"))
cbind(c(min(out_gibbs$B[1,1, ]), min(out_gibbs$B[2,1, ])),
      c(max(out_gibbs$B[1,1, ]), max(out_gibbs$B[2,1, ])))
# out_gibbs <- new_GibbsOutputIIG(out_gibbs)
summary_out <- summary(out_gibbs, true_vals = true_vals,
                       options = list(nrows_out = 9999,
                                      sig_digits = 4,
                                      roundoff = 4,
                                      ci_val = 1.96))

path_to_save_all <- "./data/output/tmp"
f_pm <- apply(out_gibbs$f[,, burnin:itermax], c(1, 2), mean)
B_pm <- apply(out_gibbs$B[,, burnin:itermax], c(1, 2), mean)
D_pm <- apply(out_gibbs$D[,, burnin:itermax], c(1, 2), mean)
plot_fit_diagnostics(yObs = yObs, wRegs = w_regs, 
                     f_pm = f_pm, B_pm = B_pm, 
                     D_pm = D_pm, NN_num_y = NN * num_y,
                     settings = list(path_save = path_to_save_all,
                                     plot_grid = c(3, 3)))
plot_factor_diagnostics(f_means_simul, f_pm,
                        NN_num_y = NN * num_y, num_joint_fac = njfacs,
                        settings = list(path_save = path_to_save_all,
                                     plot_grid = c(3, 3)))
plot_fBido_diagnostics(f_true = f_means_simul,
                       B_true = B_means,
                       f_pm = f_pm, B_pm = B_pm, 
                       NN_num_y = NN * num_y,
                       settings = list(path_save = path_to_save_all,
                                       plot_grid = c(3, 3)))
plot_fBjnt_diagnostics(f_true = f_means_simul,
                       B_true = B_means,
                       f_pm = f_pm, B_pm = B_pm, 
                       NN_num_y = NN * num_y, num_joint_fac = njfacs,
                       settings = list(path_save = path_to_save_all,
                                       plot_grid = c(3, 3)))

IIGpkg::save_Gibbs_plots(
  Gibbs = out_gibbs,
  burnin = burnin,
  path = input_simul_study$plotPath,
  nameMat = input_simul_study$nameMat,
  nameReg = input_simul_study$nameRegList[[out_gibbs$initials$wRegSpec]],
  njointfac = out_gibbs$initials$njointfac,
  type = input_simul_study$type,
  predictionCI = TRUE,
  onlyY = FALSE,
  hierachPrior = TRUE)
