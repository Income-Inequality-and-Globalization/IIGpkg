
source("GibbsSSM_2.R")







Gibbs2_SM_SA_sampler <- function(p_joint, B_par, D_par, mu_b0_par, Sigma_b0_par, alpha_b0_par, beta_b0_par, nreg, OmegaLoad0Scale, OmegaReg0Scale, countryA,  A_diag, nu0, Psi0, shape0, rate0, yObs, wRegSpec, wReg, Vhat, incObsOld = 100000, incObsNew = 100000, covScale = 1, VhatDiagScale, VhatDiagScale_start, VdiagEst, alpha0, beta0, N, TT,storePath,itermax, scaleA, diagA, sampleA, identification, type = "allidio"){

  npara <- 3
  if(type == "allidio"){
    p <- npara * N + p_joint
    # Selection matrices
    # selectR_vec <- c(rep(1, npara*p_joint+1), rep(c(rep(0,npara),1),npara - 1)) 
    selectR_vec <- c(as.numeric(c(lower.tri(matrix(rep(1,npara * p_joint), npara, p_joint),diag=T))),1,rep(c(rep(0,npara),1),npara - 1))
    selectR <- diag(selectR_vec)[-which(selectR_vec == 0),]
    selectC <- rep(0, npara + npara * p_joint)
    
    OmegaLoad0 <- OmegaLoad0Scale * diag((p_joint + npara))
    
  }else if(type == "countryidio"){
    p <- N + p_joint
    # Selection matrices
    selectR <- diag(npara * (p_joint + 1))
    selectC <- rep(0, npara * (p_joint + 1))
    
    OmegaLoad0 <- OmegaLoad0Scale * diag((p_joint + 1))
  }else if(type == "countryidio_nomu"){
    p <- N + p_joint
    # Selection matrices
    selectR <- diag(npara * (p_joint + 1))[-npara*(p_joint +1),]
    selectC <- rep(0, npara * (p_joint + 1))[-npara*(p_joint +1)]
    
    OmegaLoad0 <- OmegaLoad0Scale * diag((p_joint + 1))
  }
  
  OmegaReg0 <- OmegaReg0Scale * diag(npara * nreg)
  Omega0 <- as.matrix(Matrix::bdiag(OmegaLoad0, OmegaReg0))
  const <- rep(0, npara * N)
  #B_par <- B_par
  #B_par <- apply(yObs,1,\(x) sd(diff(x[!is.na(x)]), na.rm = T) )
  B_start <- makeBstart(npara = npara, N = N, p = p, p_joint = p_joint, B_par = B_par, type = type)
  B_start_stack <- B_start$B_stack
  B_prior_i <- B_start$B_i
  A <- A_diag*diag(npara)
  
  
  if(nreg == 0){
    wReg <- t(as.matrix(rep(1, TT)))
    D_start_stack <- NULL
    D_prior_i <- NULL
  }else{
    lr <- diag(nreg * npara)
    ur <- matrix(0, npara*(p_joint +1),npara*nreg)
    ll <- matrix(0,npara*nreg,npara^2 + npara*p_joint)
    selectR <- cbind(rbind(selectR, ll), rbind(ur,lr) )
    selectC <- rep(0, npara + npara * p_joint + npara * nreg)
    Dstart <- makeDstart(npara = npara, N = N, nreg = nreg, D_par = D_par)
    D_start_stack <- Dstart$Dstack
    D_prior_i <- Dstart$D_i
  }
  
  initf <- rep(0,p)
  initP <- diag(0,p)
  #VDiag_start <- rep(1,npara * N)

  Q <- diag(p)
  Phi <- diag(p)
  if(missing(VhatDiagScale_start)){
    VhatDiagScale_start <- NULL
  }


  Vhat <- array(apply(Vhat, 3, covarianceScale, scale = covScale), c(npara, npara, N * TT))

  ncoef <- sum(selectR)
  mu_b0 <- rep(mu_b0_par, ncoef)
  Sigma_b0 <- diag(Sigma_b0_par, ncoef)
  alpha_b0 <- rep(alpha_b0_par, ncoef)
  beta_b0 <- rep(beta_b0_par, ncoef)
  

  
  identmax <- 1000
  itermax <- itermax
  storePath <- storePath
  storeUnit <- 10000
  
  set.seed(123)
  start <- Sys.time()
  Gibbs2_SM_SA <- GibbsSSM_2(itermax= itermax, identmax = identmax, npara = npara, nreg = nreg, njointfac = p_joint, type = type, yObs = yObs, c_= const, D = D_start_stack, D0 = D_prior_i,  B = B_start_stack ,Phi=Phi,Q = Q, initX = initf, initP = initP, initU = 0, wRegSpec = wRegSpec,wReg = wReg,
                             B0 = B_prior_i, mu_b0 = mu_b0, Sigma_b0 = Sigma_b0, alpha_b0 = alpha_b0, beta_b0 = beta_b0, Omega0 = Omega0, selectR = selectR, selectC = selectC, Vhat = Vhat, incObsOld = incObsOld, incObsNew = incObsNew, covScale = covScale, VhatDiagScale = VhatDiagScale, VhatDiagScale_start = VhatDiagScale_start, VdiagEst = VdiagEst, alpha0 =alpha0, beta0 = beta0, countryA = countryA, A = A, nu0 = nu0, Psi0 = Psi0, shape0 = shape0, rate0 = rate0,
                             storePath = storePath, storeUnit =  storeUnit, scaleA = scaleA, diagA = diagA, sampleA = sampleA, identification = identification)
  # ,fPost = t(test_data$f) for checking the single parts of the sampler
  #saveRDS(test_Gibbs,"test_Gibbs_initP_0_N2_T500_100000.rds")
  #rm(test_Gibbs)
  end <- Sys.time()
  
  return(list(time_diff = end-start, Gibbs2_SM_SA =Gibbs2_SM_SA))
}




covarianceScale <- function(scale, Var){
  if(scale == 1){
   return(Var)  
  }else if(any(is.na(Var))){
    return(Var)
  }else{
    n <- ncol(Var)
    Cor <- cov2cor(Var)
    return(diag(sqrt(diag(Var))) %*% (Cor*scale + (1 - scale) * diag(n)) %*% diag(sqrt(diag(Var))))
  }
}

