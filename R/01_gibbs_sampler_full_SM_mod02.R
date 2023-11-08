#' Wrapper over the full Gibbs sampler
#'
#' @inheritParams GibbsSSM_2
#' @param p_joint Anzahl gmeinsamer Faktoren (njointfac in FFBS und den helpers)
#' @param B_par Vektor der Prior-Erwartungswerte fuer die Ladungen auf die latenten Faktoren (Wenn B_par ein Skalar ist, wird dieser Wert fuer alle Ladungen genutzt)
#' @param D_par Vektor der Prior-Erwartungswerte fuer die partiellen Effekte der Regressoren (Wenn B_par ein Skalar ist, wird dieser Wert fuer alle Ladungen genutzt)
#' @param nreg Anzahl der Regressor-Variablen
#' @param OmegaLoad0Scale Vektor der Prior-Varianzen der Ladungen (der latenten Fakotren) (wenn OmegaLoad0Scale ein Skalar ist, wird dieser Wert fuer jede Varianz genutzt)
#' @param OmegaReg0Scale Vektor der Prior-Varianzen der partiellen Effekte (der Regressoren) (wenn OmegaReg0Scale ein Skalar ist, wird dieser Wert fuer jede Varianz genutzt)
#' @param A_diag MISSING DESCRIPTION TO ADD LATER
#' @param TT Time dimension (number of points in time)
#' @param N Cross-sectional dimension (number of countries)
# @param countryA MISSING DESCRIPTION TO ADD LATER
# @param nu0 MISSING DESCRIPTION TO ADD LATER
# @param Psi0 MISSING DESCRIPTION TO ADD LATER
# @param shape0 MISSING DESCRIPTION TO ADD LATER
# @param rate0 MISSING DESCRIPTION TO ADD LATER
# @param yObs MISSING DESCRIPTION TO ADD LATER
# @param wRegSpec MISSING DESCRIPTION TO ADD LATER
# @param wReg MISSING DESCRIPTION TO ADD LATER
# @param Vhat MISSING DESCRIPTION TO ADD LATER
# @param incObsOld MISSING DESCRIPTION TO ADD LATER
# @param incObsNew MISSING DESCRIPTION TO ADD LATER
# @param covScale MISSING DESCRIPTION TO ADD LATER
# @param VhatDiagScale MISSING DESCRIPTION TO ADD LATER
# @param VhatDiagScale_start MISSING DESCRIPTION TO ADD LATER
# @param VdiagEst MISSING DESCRIPTION TO ADD LATER
# @param alpha0 MISSING DESCRIPTION TO ADD LATER
# @param beta0 MISSING DESCRIPTION TO ADD LATER
# @param storePath MISSING DESCRIPTION TO ADD LATER
# @param itermax MISSING DESCRIPTION TO ADD LATER
# @param scaleA MISSING DESCRIPTION TO ADD LATER
# @param diagA MISSING DESCRIPTION TO ADD LATER
# @param sampleA MISSING DESCRIPTION TO ADD LATER
# @param identification MISSING DESCRIPTION TO ADD LATER
# @param type MISSING DESCRIPTION TO ADD LATER
#'
#' @return Gibbs sampler output
#' @export
Gibbs2_SM_SA_sampler <- function(p_joint,
                                 B_par,
                                 D_par,
                                 prior_list = NULL,
                                 nreg,
                                 OmegaLoad0Scale,
                                 OmegaReg0Scale,
                                 countryA,
                                 A_diag,
                                 nu0,
                                 Psi0,
                                 shape0, rate0,
                                 yObs,
                                 wRegSpec,
                                 wReg,
                                 Vhat,
                                 incObsOld = 100000,
                                 incObsNew = 100000,
                                 covScale = 1,
                                 VhatDiagScale,
                                 VhatDiagScale_start,
                                 VdiagEst,
                                 alpha0,
                                 beta0,
                                 N,
                                 TT,
                                 storePath,
                                 itermax,
                                 scaleA,
                                 diagA,
                                 sampleA,
                                 identification,
                                 type = "allidio") {
  mu_b0_par    <- prior_list$hyperpriors$mu_b0_par
  Sigma_b0_par <- prior_list$hyperpriors$Sigma_b0_par
  alpha_b0_par <- prior_list$hyperpriors$alpha_b0_par
  beta_b0_par  <- prior_list$hyperpriors$beta_b0_par
  
  ## Number of parameters from GMM estimation: a, q, mu
  npara <- 3

  if (type == "allidio") {
    p <- npara * N + p_joint
    # Selection matrices
    # selectR_vec <- c(rep(1, npara*p_joint+1), rep(c(rep(0,npara),1),npara - 1))
    selectR_vec <- c(as.numeric(c(lower.tri(matrix(rep(1, npara * p_joint), npara, p_joint), diag = T))), 1, rep(c(rep(0, npara), 1), npara - 1))
    selectR <- diag(selectR_vec)[-which(selectR_vec == 0), ]
    selectC <- rep(0, npara + npara * p_joint)

    OmegaLoad0 <- OmegaLoad0Scale * diag(npara * (p_joint + npara))
  } else if (type == "countryidio") {
    p <- N + p_joint
    # Selection matrices
    selectR <- diag(npara * (p_joint + 1))
    selectC <- rep(0, npara * (p_joint + 1))

    OmegaLoad0 <- OmegaLoad0Scale * diag(npara * (p_joint + 1))
  } else if (type == "countryidio_nomu") {
    p <- N + p_joint
    # Selection matrices
    selectR <- diag(npara * (p_joint + 1))[-npara * (p_joint + 1), ]
    selectC <- rep(0, npara * (p_joint + 1))[-npara * (p_joint + 1)]

    OmegaLoad0 <- OmegaLoad0Scale * diag(npara * (p_joint + 1))
  }

  OmegaReg0 <- OmegaReg0Scale * diag(npara * nreg)
  Omega0 <- as.matrix(Matrix::bdiag(OmegaLoad0, OmegaReg0))
  const <- rep(0, npara * N)
  # B_par <- B_par
  # B_par <- apply(yObs, 1, \(x) sd(diff(x[!is.na(x)]), na.rm = T))
  B_start <- makeBstart(npara = npara, N = N, p = p, p_joint = p_joint, B_par = B_par, type = type)
  B_start_stack <- B_start$B_stack # weiterbenutzen
  B_prior_i <- B_start$B_i # weiterbenutzen
  A <- A_diag * diag(npara) # weiterbenutzen

  if (nreg == 0) {
    wReg <- t(as.matrix(rep(1, TT))) # regressoren festsetzen
    D_start_stack <- NULL
    D_prior_i <- NULL
  } else {
    lr <- diag(nreg * npara)
    ur <- matrix(0, npara * (p_joint + 1), npara * nreg)
    ll <- matrix(0, npara * nreg, npara^2 + npara * p_joint)
    selectR <- cbind(rbind(selectR, ll), rbind(ur, lr))
    selectC <- rep(0, npara + npara * p_joint + npara * nreg)
    Dstart <- makeDstart(npara = npara, N = N, nreg = nreg, D_par = D_par)
    D_start_stack <- Dstart$Dstack
    D_prior_i <- Dstart$D_i
  }

  initf <- rep(0, p)
  initP <- diag(0, p)
  initU <- 0
  # VDiag_start <- rep(1,npara * N)

  Q <- diag(p)
  Phi <- diag(p)
  if (missing(VhatDiagScale_start)) {
    VhatDiagScale_start <- NULL
  }

  Vhat <- array(apply(Vhat, 3, covarianceScale, scale = covScale), c(npara, npara, N * TT))

  identmax <- 1000
  itermax <- itermax
  storePath <- storePath
  storeUnit <- 10000

  prior_list$priors$B0 <- B_prior_i
  prior_list$priors$Omega0 <- Omega0
  start <- Sys.time()
  set.seed(123)
  Gibbs2_SM_SA <- GibbsSSM_2(
    itermax = itermax,
    identmax = identmax,
    npara = npara,
    nreg = nreg,
    njointfac = p_joint,
    type = type,
    yObs = yObs,
    c_ = const,
    D = D_start_stack,
    D0 = D_prior_i,
    B = B_start_stack,
    Phi = Phi,
    Q = Q,
    initX = initf,
    initP = initP,
    initU = initU,
    wRegSpec = wRegSpec,
    wReg = wReg,
    prior_list = prior_list,
    selectR = selectR,
    selectC = selectC,
    Vhat = Vhat,
    incObsOld = incObsOld,
    incObsNew = incObsNew,
    covScale = covScale,
    VhatDiagScale = VhatDiagScale,
    VhatDiagScale_start = VhatDiagScale_start,
    VdiagEst = VdiagEst,
    alpha0 = alpha0,
    beta0 = beta0,
    countryA = countryA,
    A = A,
    nu0 = nu0,
    Psi0 = Psi0,
    shape0 = shape0,
    rate0 = rate0,
    storePath = storePath,
    storeUnit = storeUnit,
    scaleA = scaleA,
    diagA = diagA,
    sampleA = sampleA,
    identification = identification
  )
  # ,fPost = t(test_data$f) for checking the single parts of the sampler
  # saveRDS(test_Gibbs,"test_Gibbs_initP_0_N2_T500_100000.rds")
  # rm(test_Gibbs)
  end <- Sys.time()

  return(list(time_diff = end - start, Gibbs2_SM_SA = Gibbs2_SM_SA))
}



#' Covariance/Corrleation scale
#'
#' @param scale Scaling factor
#' @param Var VCOV matrix
#'
#' @return returns covariance/correlation scale
#' @export
covarianceScale <- function(scale, Var) {
  if (scale == 1) {
    return(Var)
  } else if (any(is.na(Var))) {
    return(Var)
  } else {
    n <- ncol(Var)
    Cor <- cov2cor(Var)
    return(diag(sqrt(diag(Var))) %*% (Cor * scale + (1 - scale) * diag(n)) %*% diag(sqrt(diag(Var))))
  }
}
