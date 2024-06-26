#' MCMC sampler for hyperprior-mean
#'
#' @param NN cross sectional units
#' @param npara dimension of measurements (of income distribution)
#' @param nreg number of regressors
#' @param njointfac number of joint factors
#' @param B_i B-matrices factorized per i
#' @param D_i D-matrices factorized per i
#' @param invOmega0 inverse of prior matrix
#' @param mu_b0 hyperprior parameter
#' @param Sigma_b0 hyperprior parameter
#' @param selectR selection matrix
#' @param type sampling type
#'
#' @return a single posterior sample for hyperprior-parameter mu
#' @export
mu_sampler <- function(NN, npara, nreg,
                       njointfac,
                       B_i, D_i,
                       invOmega0, mu_b0, Sigma_b0,
                       selectR, type) {
  BD <- sapply(1:NN, \(x) selectR %*% c(B_i[, , x], D_i[, , x]))
  meanBD <- apply(BD, 1, mean)
  invSigma_b0 <- solve(Sigma_b0)
  Sigma_b1_inv <- invSigma_b0 + NN * invOmega0
  Sigma_b1 <- solve((Sigma_b1_inv + t(Sigma_b1_inv)) / 2)
  mu_b1 <- Sigma_b1 %*% (invSigma_b0 %*% mu_b0 + NN * invOmega0 %*% meanBD)

  mu0_sim <- mvtnorm::rmvnorm(n = 1, mean = mu_b1, sigma = Sigma_b1)

  B0vec <- mu0_sim[1:((njointfac + 1) * npara)]

  if(nreg == 0){
    D0 <- NULL
  } else {
    D0vec <- mu0_sim[-(1:((njointfac + 1) * npara))]
    D0 <- replicate(NN, matrix(D0vec, nrow = npara, ncol = nreg))
  }


  if (njointfac == 1) {
    idioFac <- B0vec[-(1:npara)]
    # if (type == "countryidio_nomu") {
    #   idioFac <- c(idioFac, 0)
    # }
    jointFac <- matrix(B0vec[1:(npara * njointfac)], ncol = njointfac)
    B0 <- cbind(jointFac, diag(idioFac))
  } else if (njointfac == 0) {
    idioFac <- B0vec
    B0 <- diag(idioFac)
  } else {
    stop("Not yet implemented.")
  }

  B0 <- replicate(NN, B0)
  return(list(B0 = B0, D0 = D0))
}
#' MCMC sampler for hyperprior-VCM
#'
#' @inheritParams mu_sampler
#' @param alpha_b0 hyperprior parameter
#' @param beta_b0 hyperprior parameter
#'
#' @return a single posterior sample for hyperprior-parameter Omega
#' @export
omega_sampler <- function(NN, B_i, D_i, B0, D0, alpha_b0, beta_b0, selectR) {
  BD <- sapply(1:NN, \(x) selectR %*% c(B_i[, , x], D_i[, , x]))
  BD_0 <- c(B0[, , 1], D0[, , 1])

  k <- sum(selectR)
  Omega0 <- numeric(k)
  for (i in 1:k) {
    alpha_b1 <- alpha_b0[i] + NN / 2
    beta_b1 <- beta_b0[i] + sum((BD[i] - BD_0[i])^2) / 2

    Omega0[i] <- LaplacesDemon::rinvgamma(1, shape = alpha_b1, scale = beta_b1)
  }

  diag(Omega0)
}
