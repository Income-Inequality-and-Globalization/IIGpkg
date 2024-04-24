rm(list = ls())
library(tidyverse)
library(IIGpkg)
# source("plot_functions_ACF.R")
# Data --------------------------------------------------------------------
# GMM Ergebnisse (Parameter-Schaetzungen) mit Jahr 2021
pth_base_data    <- "./data/input/data-sm"
# data_raw_update_and_save(pth_base_data,
#                      data_file = "data_coef_covariates.txt",
#                      save_file = "data_coef_covariates_2020.rds",
#                      year_to_remove = 2021)
# GMM Ergebnisse (Kovarianzen) nur bis 2020
# pth_data_covr           <- file.path(pth_base_data, "data_coef_covariates_2020.rds")
pth_data_covr           <- file.path(pth_base_data, "data_coef_covariates_2020_grw.rds")
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

# regs <- c("cpi_change", "unemployment", "gdp_ppp")
regs <- c("cpi_change", "unemployment", "gdp_ppp_pc")
# regs <- c("cpi_change", "unemployment", "gdp_ppp_2017_grw")
# GMM_by_year$gdp_ppp_pc <- GMM_by_year$gdp_ppp / GMM_by_year$population
# GMM_by_year$gdp_ppp_pc[GMM_by_year$country == "Venezuela"] <- GMM_by_year$gdp_ppp[GMM_by_year$country == "Venezuela"]
tmplist_rg <- generate_regressor_combinations(GMM_by_year,
                                              regs,
                                              TT,
                                              countries,
                                              years)
nameRegList   <- tmplist_rg$nameRegList
wRegRawList   <- tmplist_rg$wRegRawList
wRegCombsList <- tmplist_rg$wRegCombsList


wRegCombsList[[8]] <- matrix(0,0,0)






shape0 <- 2
rate0 <- 1
type <- "allidio"
alpha0 <- 0.5
beta0 <- 0.5




itermax <- 200
burnin  <- 100


path <- "./data/output/empirical"
path <- "D:/IIGpkg_output_test"

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


############## Cluster ##############
p_joint_val <- c(0,1)
reg_spec_val <- c(7, 8)

n_cluster <- 4
grid_ <- expand.grid(p_joint_val, reg_spec_val)
index_grid <- matrix(1:dim(grid_)[1], ncol = n_cluster, byrow=T)

#####################################



for(j in 1:dim(index_grid)[1]){

  cl <- parallel::makeCluster(n_cluster)

  parallel::clusterExport(cl, c("grid_", "index_grid", "spec_01", "shape0", "rate0","yObs_log_centered_standardized", "wRegCombsList", "VdiagEst", "VhatDiagScale", "VCOV_array_country_pd", "alpha0", "beta0", "N", "TT",
                                "storePath", "itermax", "type"))



  #out <-  parallel::parApply(cl, grid_[index_grid[j,],], 1, \(x) as.numeric(x[2]))


  out <-  parallel::parApply(cl, grid_[index_grid[j,],], 1, \(x){IIGpkg::Gibbs2_SM_SA_sampler(
               p_joint = as.numeric(x[1]),
                B_par = 0,
                D_par = 0, #grid[j, 6],
                prior_list = spec_01$prior_list,
                nreg = ifelse(length(dim(wRegCombsList[[as.numeric(x[2])]])[1]/N) != 0, dim(wRegCombsList[[as.numeric(x[2])]])[1]/N, 0  ),
                OmegaLoad0Scale = spec_01$OmegaLoad0Scale, #grid[jj, 1],
                OmegaReg0Scale = spec_01$OmegaReg0Scale, #grid[jj, 2],
                countryA = FALSE,
                A_diag = 1,
                nu0 = spec_01$nu0, # 7, #grid[jj, 3],  #7, #grid[jj, 3],
                Psi0 = spec_01$Psi0, # diag(3) * 57, # grid[jj, 4], # 27,# grid[jj, 4],
                shape0 = shape0,
                rate0 = rate0,
                yObs = yObs_log_centered_standardized,
                wRegSpec = as.numeric(x[2]),
                wReg = wRegCombsList[[as.numeric(x[2])]],
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
                storePath = storePath, #".", #storePath,
                itermax = itermax,
                scaleA = FALSE,
                diagA = FALSE,
                sampleA = TRUE,
                identification = TRUE,
                type = type)}
                )


  parallel::stopCluster(cl)

}



# out <- IIGpkg::Gibbs2_SM_SA_sampler(
#   p_joint = 0, #p_joint,
#   B_par = 0,
#   D_par = 0, #grid[j, 6],
#   prior_list = spec_01$prior_list,
#   nreg = 0,#dim(wRegCombsList[[7]])[1]/N,
#   OmegaLoad0Scale = spec_01$OmegaLoad0Scale, #grid[jj, 1],
#   OmegaReg0Scale = spec_01$OmegaReg0Scale, #grid[jj, 2],
#   countryA = FALSE,
#   A_diag = 1,
#   nu0 = spec_01$nu0, # 7, #grid[jj, 3],  #7, #grid[jj, 3],
#   Psi0 = spec_01$Psi0, # diag(3) * 57, # grid[jj, 4], # 27,# grid[jj, 4],
#   shape0 = shape0,
#   rate0 = rate0,
#   yObs = yObs_log_centered_standardized,
#   wRegSpec = 7,
#   wReg = wRegCombsList[[7]],
#   Vhat = VCOV_array_country_pd,
#   incObsOld = 100000,
#   incObsNew = 5000,
#   covScale = 1,
#   VdiagEst = VdiagEst,
#   VhatDiagScale = VhatDiagScale,
#   alpha0 = alpha0,
#   beta0 = beta0,
#   N = N,
#   TT = TT,
#   storePath = storePath, #".", #storePath,
#   itermax = itermax,
#   scaleA = FALSE,
#   diagA = FALSE,
#   sampleA = TRUE,
#   identification = TRUE,
#   type = type)
