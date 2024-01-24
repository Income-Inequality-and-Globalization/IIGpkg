library(tidyverse)
library(IIGpkg)
# source("plot_functions_ACF.R")
# Data --------------------------------------------------------------------
firstObs_center <- function(x){
  nonNA_index <- which(!is.na(x))[1]
  x - x[nonNA_index]
}
firstObs_center_values <- function(x){
  x[which(!is.na(x))[1]]
}
standardize <- function(x){
  x/sd(x, na.rm = TRUE)
}
standardized_values <- function(x){
  sd(x, na.rm = TRUE)
}


# GMM Ergebnisse (Parameter-Schaetzungen) mit Jahr 2021
pth_base_data <- "./data/input/data-sm"
pth_data_covr <- file.path(pth_base_data, "data_coef_covariates.txt")
results_GMM_2021 <- tibble(read.table(pth_data_covr,header = TRUE))
# GMM Ergebnisse (Kovarianzen) mit Jahr 2021
VCOV_array_country_2021 <- readRDS("./data/input/VarianceArray_by_country.rds")
# 2021 wird entfernt, das es fuer alle Laender fehlt.
#### aus Datensatz entefernen und Befehl loeschen
VCOV_array_country <- VCOV_array_country_2021[, , which(results_GMM_2021$year != 2021)] 

GMM_by_year <- results_GMM_2021 %>% arrange(year) %>% filter(year != 2021)
years <- unique(GMM_by_year$year)
TT <- length(years)
npara <- 3

countries <- unique(GMM_by_year$country)
N <- length(countries)
nameMat <- cbind(rep(countries, each = npara), rep(c("a","q","mu"), N))

VCOV_array_country_pd <- array(apply(VCOV_array_country,3,function(x){x/2 + t(x)/2 }),  c(npara, npara, N * TT)) 

Diag_VCOV <- FALSE

if (Diag_VCOV) {
  VCOV_array_country_pd <-  array(apply(VCOV_array_country_pd,3,\(x) diag(diag(x)) ), c(npara, npara, N * TT))
}
yObs <- matrix(t(GMM_by_year %>%  select(a,q,est_mean)), ncol= TT)
rownames(yObs) <- rep(countries, each = npara)
colnames(yObs) <- years

regs <- c("cpi_change", "unemployment", "gdp_ppp")
lregs <- length(regs) # Anzahl der verfuegbaren Regressoren (= 3)
reg_combs <- lapply(1:lregs, \(x) combn(regs, x)) # Alle Kombinationsmoeglichkeiten der Regressoren
ncombs <- sum(sapply(1:lregs, \(x)choose(3,x))) # Anzahl der Regressorkombinationen (= 7)
wRegCombsList <- vector("list", ncombs) # Liste der Regressormatrizen fuer alle Kombinationen (wird in der foor-loop befüllt)
nameRegList <- vector("list", ncombs) # Alle Kombinationsmoeglichkeiten der Regressoren (wie reg_combs, aber als Liste mit 7 Elementen) (wird in der foor-loop befüllt)

ii <- 1
for (i in 1:length(reg_combs)){
   for (j in 1:ncol(reg_combs[[i]])){
      wReg <- matrix(t(GMM_by_year %>%  select(paste(reg_combs[[i]][,j]) )), ncol= TT)
      wReg_centered <- t(apply(wReg,1,firstObs_center))
      wReg <- t(apply(wReg_centered,1,standardize))
      wReg[is.na(wReg)] <- 0
      rownames(wReg) <- rep(countries, each = nrow(reg_combs[[i]]))
      colnames(wReg) <- years
      nameRegList[[ii]] <- as.character(sapply(paste(reg_combs[[i]][,j]), \(x) substr(x,1,7) ))
      wRegCombsList[[ii]] <- wReg
      ii <- ii + 1
      rm(wReg)
   }
}
rm(ii)





  

# Log ---------------------------------------------------------------------

yObs_unadj <- yObs
yObs <- log(yObs)

VCOV_logAdj <- function(npara, N, TT, yObs, VCOV_array){
  yObs_inv <- 1/yObs
  ArrayTimeCountry <- array(yObs_inv, c(npara,1, N*TT ))
  MatCountryTime <-  ArrayTimeCountry[,,c(sapply(0:(N-1),\(x)(x+seq(1,(TT-1)*N+1,N))))]
  VCOV_log_list <- lapply(1:(N*TT),\(x) diag(MatCountryTime[,x]) %*% VCOV_array[,,x] %*% diag(MatCountryTime[,x])  ) 
  VCOV_log_array <-  array(unlist(VCOV_log_list), dim = c(npara,npara,TT*N))
  array(apply(VCOV_log_array,3,function(x){x/2 + t(x)/2 }),  c(npara, npara, N * TT))
}

VCOV_array_country_pd <- VCOV_logAdj(npara, N, TT, yObs_unadj, VCOV_array_country_pd)



sd_yObs <- apply(yObs, 1, sd, na.rm = TRUE)

standardize_VCOV <- function(VCOV_array, TT, N, sd_vec, npara){
  V_adj_array <- array(0, dim = c(npara, npara, TT*N))
  
  for(i in 1:N){
    V_adj <- diag(1 / sd_vec[((i-1)*npara + 1):(i * npara )])
    V_select <- VCOV_array[,,((i-1) * TT + 1):(i * TT)]
    adj_array <-  array(apply(V_select, 3, \(x) V_adj %*% x %*% V_adj), dim = c(npara, npara, TT))
    V_adj_array[,,((i-1) * TT + 1):(i * TT)] <- adj_array
  }
  
  return(V_adj_array)
}

VCOV_array_country_pd <- standardize_VCOV(VCOV_array_country_pd, TT, N, sd_yObs, npara)


yObs_centered <- t(apply(yObs, 1, firstObs_center))
yObs_centered_values <- t(apply(yObs, 1, firstObs_center_values))
yObs <- t(apply(yObs_centered, 1, standardize))
yObs_standardization <- t(apply(yObs_centered, 1, standardized_values))



n_cluster <- 12
#B_par_values <- c(1, 2, 5, 10)
#B_par_values <- c(2)
B_par_values <- c(0.1,1,5)
Omega0LoadScale_values <- c(0.001,0.1,1,5)

D_par_values <- c(0, 1, 5)
Omega0RegScale_values <- c(0.001,1,5)

sampleA_values <- c(F,T)
#p_joint <- c(0,1)
#covScale_values <- c(0, 0.01, 0.2,0.5,1)
covScale_values <- c(1,25,50,100)
Adiag_values <- c(F,T)

wReg_values <- c(7)

nu0 <- 7
nu0_values <- c(7, 60, 500, 1000)
Psi0 <- 2.5*diag(3)
Psi0_values <- c(1,500,1000) 
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

#grid <- expand.grid(B_par_values, Omega0Scale_values, A_diag_values)
#grid <- expand.grid(B_par_values, Omega0Scale_values, sampleA_values)
#grid <- expand.grid(B_par_values, Omega0Scale_values, covScale_values)
B_par_values <- 1
D_par_values <- 1
wReg_values <- 7
Omega0LoadScale_values <- 1000
Omega0RegScale_values <- 1000
grid <- expand.grid(Omega0LoadScale_values, Omega0RegScale_values, nu0_values,Psi0_values, wReg_values, D_par_values,B_par_values)
index_grid <- matrix(1:dim(grid)[1], ncol=n_cluster,byrow=T)


itermax <- 200000
burnin <- 100000

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
if(VdiagEst & !VhatDiagScale){VCOV_array_country_pd <- array(rep(diag(npara), N * TT ), c(npara, npara, N*TT))}

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
#   c("yObs", "wRegCombsList", "grid", "jj", "shape0", "rate0", 
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
#       yObs = yObs,
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
      yObs = yObs,
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
      storePath = storePath,
      itermax = itermax,
      scaleA = FALSE,
      diagA = FALSE,
      sampleA = TRUE,
      identification = TRUE,
      type = type)
IIGpkg::save_Gibbs_plots(
  Gibbs = out$Gibbs2_SM_SA,
  burnin = burnin,
  path = plotPath,
  nameMat = nameMat,
  nameReg = nameRegList[[out$Gibbs2_SM_SA$initials$wRegSpec]],
  njointfac = out$Gibbs2_SM_SA$initials$njointfac,
  type = type,
  predictionCI = TRUE,
  onlyY = FALSE,
  hierachPrior = TRUE)
