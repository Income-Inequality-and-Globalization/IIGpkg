rm(list = ls())
library(tidyverse)
library(IIGpkg)
# source("plot_functions_ACF.R")
# Data --------------------------------------------------------------------
# GMM Ergebnisse (Parameter-Schaetzungen) mit Jahr 2021
pth_base_data    <- "./data/input/data-sm"
# remove_year_and_save(pth_base_data,
#                      data_file = "data_coef_covariates.txt",
#                      save_file = "data_coef_covariates_2020.rds",
#                      year_to_remove = 2021)
# GMM Ergebnisse (Kovarianzen) nur bis 2020
pth_data_covr           <- file.path(pth_base_data, "data_coef_covariates_2020.rds")
pth_data_vcov           <- file.path(pth_base_data, "VarianceArray_by_country_2020.rds")
GMM_by_year             <- readRDS(file = pth_data_covr)
VCOV_array_country_2020 <- readRDS(pth_data_vcov)

data_meta_info <- get_data_meta_attributes(GMM_by_year)
countries      <- data_meta_info$countries
nameMat        <- data_meta_info$nameMat
npara          <- data_meta_info$npara
years          <- data_meta_info$years
TT             <- data_meta_info$TT
N              <- data_meta_info$N

VCOV_array_country_pd <- process_covariance_array(VCOV_array_country_2020,
                                                  npara,
                                                  N,
                                                  TT)
tmplist_dt <- process_observation_data(GMM_by_year, 
                                       npara,
                                       TT,
                                       countries,
                                       years)
# Get various data sets: un-adjusted version, in logs, in logs and centered, 
# and finally in logs, centered and standardized:
yObs_unadj                     <- tmplist_dt$data$yObs_unadj
yObs_log                       <- tmplist_dt$data$yObs_log
yObs_log_centered              <- tmplist_dt$data$yObs_log_centered
yObs_log_centered_standardized <- tmplist_dt$data$yObs_log_centered_standardized
# Get the values used for above transformations
yObs_log_centered_values <- tmplist_dt$data_values$yObs_log_centered_values
yObs_standardization     <- tmplist_dt$data_values$yObs_standardized_values 
sd_yObs_log              <- tmplist_dt$data_values$sd_yObs_log
# Adjust the measurement VCOV matrix (`A`) to this log- and other transformation
VCOV_array_country_pd <- VCOV_array_country_pd %>%
  log_adj_VCOV(npara, N, TT, yObs_unadj) %>%
  standardize_VCOV(TT, N, sd_yObs_log, npara)

regs <- c("cpi_change", "unemployment", "gdp_ppp")
tmplist_rg <- generate_regressor_combinations(GMM_by_year,
                                              regs,
                                              TT,
                                              countries,
                                              years) 
nameRegList   <- tmplist_rg$nameRegList
wRegCombsList <- tmplist_rg$wRegCombsList



n_cluster <- 12
#B_par_values <- c(1, 2, 5, 10)
#B_par_values <- c(2)
B_par_values <- c(0.1, 1, 5)
Omega0LoadScale_values <- c(0.001, 0.1, 1, 5)

D_par_values <- c(0, 1, 5)
Omega0RegScale_values <- c(0.001, 1, 5)

sampleA_values <- c(FALSE, TRUE)
#p_joint <- c(0,1)
#covScale_values <- c(0, 0.01, 0.2,0.5,1)
covScale_values <- c(1, 25, 50, 100)
Adiag_values <- c(FALSE, TRUE)

wReg_values <- c(7)

nu0 <- 7
nu0_values <- c(7, 60, 500, 1000)
Psi0 <- 2.5 * diag(3)
Psi0_values <- c(1, 500, 1000) 
IncObsNew_values <- c(5000)

# mu_b0_par <- 0 
# Sigma_b0_par <- 1

# beta_b0_par <- 5
# alpha_b0_par <-  1 + beta_b0_par

#beta_b0_par/(alpha_b0_par - 1)  # IG Mean
#beta_b0_par^2/((alpha_b0_par - 1)^2*(alpha_b0_par -2)) # IG Variance


shape0 <- 2
rate0 <- 1
type <- "allidio"
alpha0 <- 0.5
beta0 <- 0.5

p_joint <- 0

# grid <- expand.grid(B_par_values, Omega0Scale_values, A_diag_values)
# grid <- expand.grid(B_par_values, Omega0Scale_values, sampleA_values)
# grid <- expand.grid(B_par_values, Omega0Scale_values, covScale_values)
B_par_values <- 1
D_par_values <- 1
wReg_values <- 7
Omega0LoadScale_values <- 1000
Omega0RegScale_values <- 1000
grid <- expand.grid(Omega0LoadScale_values,
                    Omega0RegScale_values,
                    nu0_values,
                    Psi0_values,
                    wReg_values,
                    D_par_values,
                    B_par_values)
index_grid <- matrix(1:dim(grid)[1], ncol = n_cluster,byrow = TRUE)


itermax <- 20
burnin <- 1

# storePath_vec <- c("D:/Gibbs_SM_SA_Tower/Results/Vhat_countryA_stand/Estimates", "D:/Gibbs_SM_SA_Tower/Results/Vhat_countryAdiag_stand/Estimates") 
# plotPath_vec <- c("D:/Gibbs_SM_SA_Tower/Results/Vhat_countryA_stand/Plots", "D:/Gibbs_SM_SA_Tower/Results/Vhat_countryAdiag_stand/Plots")

# path <- "D:/Gibbs_SM_SA_Tower/Results/"
path <- "./data/output/empirical"
# path <- "./data/output/simulations"
storePath <- file.path(path, "Estimates")
dir.create(storePath, recursive = TRUE)
plotPath <- file.path(path, "Plots")
dir.create(plotPath)

RegIncl <- TRUE
VdiagEst <- FALSE
VhatDiagScale <- FALSE
if (VdiagEst & !VhatDiagScale) {
  VCOV_array_country_pd <- array(rep(diag(npara), N * TT ),
                                 c(npara, npara, N*TT))
}

set.seed(123)
# cluster_specifications <- list(
#   spec_01 = list(
#     prior_list = NULL,
#     OmegaLoad0Scale = 1000,
#     OmegaReg0Scale = 1000,
#     nu0 = 7,
#     Psi0 = diag(3) * 27
#   ),
#   spec_02 = list(
#     prior_list = NULL,
#     OmegaLoad0Scale = 1000,
#     OmegaReg0Scale = 1000,
#     nu0 = 7,
#     Psi0 = diag(3) * 57
#   )
# )
# prior_list <- NULL
# n_cluster <- 2
# cl <- parallel::makeCluster(n_cluster)
jj <- 12
# parallel::clusterExport(
#   cl,
#   c("yObs_log_centered_standardized", "wRegCombsList", "grid", "jj", "shape0", "rate0", 
#     "VCOV_array_country_pd", "VhatDiagScale", "VdiagEst",
#     "alpha0", "beta0",
#     "N","TT",
#     "storePath", "plotPath", "nameMat", "nameRegList",
#     "itermax", "burnin",
#     "type",
#     "burnin")
# )
# out <- parallel::parLapply(
#   cl,
#   cluster_specifications,
#   function(x) {
#     out <- IIGpkg::Gibbs2_SM_SA_sampler(
#       p_joint = 1, #p_joint,
#       B_par = 0,
#       D_par = 0, #grid[j, 6],
#       prior_list = x$prior_list,
#       nreg = dim(wRegCombsList[[grid[jj, 5]]])[1]/N,
#       OmegaLoad0Scale = x$OmegaLoad0Scale, #grid[jj, 1],
#       OmegaReg0Scale = x$OmegaReg0Scale, #grid[jj, 2],
#       countryA = FALSE,
#       A_diag = 1,
#       nu0 = x$nu0, # 7, #grid[jj, 3],  #7, #grid[jj, 3],
#       Psi0 = x$Psi0, # diag(3) * 57, # grid[jj, 4], # 27,# grid[jj, 4],
#       shape0 = shape0,
#       rate0 = rate0,
#       yObs = yObs_log_centered_standardized,
#       wRegSpec = grid[jj, 5],
#       wReg = wRegCombsList[[grid[jj, 5]]],
#       Vhat = VCOV_array_country_pd,
#       incObsOld = 100000,
#       incObsNew = 5000,
#       covScale = 1,
#       VdiagEst = VdiagEst,
#       VhatDiagScale = VhatDiagScale,
#       alpha0 = alpha0,
#       beta0 = beta0,
#       N = N,
#       TT = TT,
#       storePath = storePath,
#       itermax = itermax,
#       scaleA = FALSE,
#       diagA = FALSE,
#       sampleA = TRUE,
#       identification = TRUE,
#       type = type)
#     IIGpkg::save_Gibbs_plots(
#       Gibbs = out$Gibbs2_SM_SA,
#       burnin = burnin,
#       path = plotPath,
#       nameMat = nameMat,
#       nameReg = nameRegList[[out$Gibbs2_SM_SA$initials$wRegSpec]],
#       njointfac = out$Gibbs2_SM_SA$initials$njointfac,
#       type = type,
#       predictionCI = TRUE,
#       onlyY = FALSE,
#       hierachPrior = TRUE)
#     return(out)
#   }
# )


spec_01 <-  list(
   prior_list = list(
      hyperpriors = list(
        mu_b0 = 0,
        Sigma_b0 = 100,
        alpha_b0 = 5, # 3, 5,  10, 255
        beta_b0 = 6 # 10, 6, 0.1, 256
      )
    ),
    OmegaLoad0Scale = 1000,
    OmegaReg0Scale = 1000,
    nu0 = 7, # 1000,
    Psi0 = diag(3) * 27 # diag(3) * 27, 57, 1000
)
out <- IIGpkg::Gibbs2_SM_SA_sampler(
      p_joint = 1, #p_joint,
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
      yObs = yObs_log_centered_standardized,
      wRegSpec = grid[jj, 5],
      wReg = wRegCombsList[[grid[jj, 5]]],
      Vhat = VCOV_array_country_pd,
      incObsOld = 100000,
      incObsNew = 5000,
      covScale = 1,
      VdiagEst = VdiagEst,
      VhatDiagScale = VhatDiagScale,
      alpha0 = alpha0,
      beta0 = beta0,
      N = N,
      TT = TT,
      storePath = ".", #storePath,
      itermax = itermax,
      scaleA = FALSE,
      diagA = FALSE,
      sampleA = TRUE,
      identification = TRUE,
      type = type)

testme <- readRDS("./cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
test1 <- identical(out$Gibbs2_SM_SA$f, testme$f)
test2 <- identical(out$Gibbs2_SM_SA$B, testme$B)
test3 <- identical(out$Gibbs2_SM_SA$D, testme$D)
test4 <- identical(out$Gibbs2_SM_SA$A, testme$A)
test5 <- identical(out$Gibbs2_SM_SA$V, testme$V)
test6 <- identical(out$Gibbs2_SM_SA$BD0STORE, testme$BD0STORE)
test7 <- identical(out$Gibbs2_SM_SA$Omega0STORE, testme$Omega0STORE)
stopifnot(all(c(test1, test2, test3, test4, test5, test6, test7)))
rm(list = c("test1", "test2", "test3", "test4", "test5", "test6", "test7"))
# IIGpkg::save_Gibbs_plots(
#   Gibbs = out$Gibbs2_SM_SA,
#   burnin = burnin,
#   path = plotPath,
#   nameMat = nameMat,
#   nameReg = nameRegList[[out$Gibbs2_SM_SA$initials$wRegSpec]],
#   njointfac = out$Gibbs2_SM_SA$initials$njointfac,
#   type = type,
#   predictionCI = TRUE,
#   onlyY = FALSE,
#   hierachPrior = TRUE)
