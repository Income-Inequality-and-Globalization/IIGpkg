#' Sum in the inverse of the posterior variance for loadings/partial effects
#'
#' @param availableObs available observations
#' @param npara number of params (usually 3)
#' @param nreg  number of regressors
#' @param njointfac number of joint factors
#' @param i cross-sectional index
#' @param fPost backward sampled states (FFBS output)
#' @param wReg regressors
#' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
#' @param type either of "allidio", "countryidio" or "countryidio_nomu"
#'
#' @return  cronecker summation part
#' @export
sumffkronV <- function(availableObs, npara, nreg, njointfac, i, fPost, wReg, Viarray, type) {
  summ <- 0
  for (tt in availableObs) {
    if (type == "allidio") {
      if (njointfac != 0) {
        f <- fPost[c(1:njointfac, (njointfac + 1 + npara * (i - 1)):(njointfac + npara * i)), tt] # generalize of no njointfac is included
      } else {
        f <- fPost[(njointfac + 1 + npara * (i - 1)):(njointfac + npara * i), tt]
      }
    } else if (type == "countryidio" | type == "countryidio_nomu") {
      if (njointfac != 0) {
        f <- fPost[c(1:njointfac, njointfac + i), tt] # generalize of no njointfac is included
      } else {
        f <- fPost[i, tt]
      }
    } else {
      stop("No valid model type selected.")
    }
    
    if (nreg != 0) {
      wReg_it <- wReg[(1 + nreg * (i - 1)):(i * nreg), tt]
      f <- c(f, wReg_it)
    }
    
    
    V <- Viarray[, , tt]
    summ <- summ + kronecker(f %*% t(f), solve(V))
  }
  return(summ)
}


#' Sum part of posterior mean for vectorized loadings/partial effects
#'
#' @inheritParams GibbsSSM_2
#' @inheritParams sumffkronV
#' @param yiObs Matrix (npara x TT) of cross-sectional unit i
#' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
#'
#' @return summation part
#' @export
sumfyV <- function(availableObs, npara, nreg, njointfac, i, fPost, wReg, yiObs, Viarray, type) {
  summ <- 0
  for (tt in availableObs) {
    if (type == "allidio") {
      if (njointfac != 0) {
        f <- fPost[c(1:njointfac, (njointfac + 1 + npara * (i - 1)):(njointfac + npara * i)), tt] # generalize of no njointfac is included
      } else {
        f <- fPost[(njointfac + 1 + npara * (i - 1)):(njointfac + npara * i), tt]
      }
    } else if (type == "countryidio" | type == "countryidio_nomu") {
      if (njointfac != 0) {
        f <- fPost[c(1:njointfac, njointfac + i), tt] # generalize of no njointfac is included
      } else {
        f <- fPost[i, tt]
      }
    } else {
      stop("No valid model type selected.")
    }
    
    if (nreg != 0) {
      wReg_it <- wReg[(1 + nreg * (i - 1)):(i * nreg), tt]
      f <- c(f, wReg_it)
    }
    
    y <- yiObs[, tt]
    V <- Viarray[, , tt]
    summ <- summ + solve(V) %*% y %*% t(f)
  }
  return(summ)
}

#' Matrix square root (by spectral decomposition)
#'
#' @param V Normal matrix
#'
#' @return matrix square root of V
#' @export
matSqrt <- function(V) {
  npara <- dim(V)[1]
  if (any(is.na(V))) {
    matrix(NA, nrow = npara, ncol = npara)
  } else {
    eig <- eigen(V)
    val <- eig$values
    vec <- eig$vectors
    # vec %*% diag(sqrt(val)) %*% t(vec)
    vec %*% diag(sqrt(val))
  }
}
#' Sorts VCOV array by time
#'
#' @inheritParams GibbsSSM_2
#' @inheritParams sumffkronV
#' @param VhatArray_A VCOV array sorted by cross-section.
#' @param N cross sectional dimension
#' @param TT time series length
#' @param Nnpara N * npara
#'
#' @return Sorted VCOV array by time
#' @export
bdiagByTime <- function(VhatArray_A, npara, N, TT, Nnpara) {
  VhatList <- plyr::alply(VhatArray_A, 3)
  VhatList_split <- split(VhatList, rep(1:TT, N))
  VhatList_bdiag <- lapply(VhatList_split, Matrix::bdiag)
  VhatArray_bdiag <- array(unlist(lapply((VhatList_bdiag), as.matrix)), c(Nnpara, Nnpara, TT))
  
  return(VhatArray_bdiag)
}
#' Sum of squared residuals (u * u') for VCOV/adjustment matrix sampling 
#'
#' @inheritParams GibbsSSM_2
#' @inheritParams bdiagByTime
#' @param uSplit List (N elements) of resiudal matrices (npara x TT) for each cross-sectional unit
#' @param VhatSqrt Array (npara x npara x N*TT) of square-root VCOVs (sorted by cross-sectional unit)
#'
#' @return sum of squared residuals
#' @export
utuSum <- function(uSplit, VhatSqrt, TT, N, npara) {
  sumUtu <- 0
  utildeSplit <- lapply(1:N, matrix, data = NA, nrow = npara, ncol = TT)
  sumUtu_individ <- lapply(1:N, matrix, data = NA, nrow = npara, ncol = npara)
  for (i in 1:N) {
    for (t in 1:TT) {
      utildeSplit[[i]][, t] <- solve(VhatSqrt[, , (i - 1) * TT + t]) %*% uSplit[[i]][, t]
    }
    
    utu <- apply(utildeSplit[[i]], 2, function(x) {
      x %*% t(x)
    })
    sumUtu_i <- apply(utu, 1, sum, na.rm = TRUE)
    sumUtu <- sumUtu + sumUtu_i
    sumUtu_i_mat <- matrix(sumUtu_i, ncol = npara)
    sumUtu_individ[[i]] <- (t(sumUtu_i_mat) + sumUtu_i_mat) / 2
  }
  return(list(sumUtu_total = matrix(sumUtu, ncol = npara), sumUtu_individ = sumUtu_individ))
}
#' Loading matrix and array
#'
#' @inheritParams Gibbs2_SM_SA_sampler
#' @inheritParams sumffkronV
#' @inheritParams bdiagByTime
#' @param p describe later
#' @param p_joint describe later 
#' @param B_par describe later
#'
#' @return Loading matrix and B for start
#' @export
makeBstart <- function(npara, N, p, p_joint, B_par, type = "allidio") { # type = c("allidio","countryidio")
  
  # if(length(B_par) == 1){B_par <- rep(B_par,N*npara*(p_joint+1))}
  B_stack <- matrix(rep(0, N * npara * p), ncol = p)
  
  if (type == "allidio") {
    non_zeros <- cbind(matrix(TRUE, ncol = p_joint, nrow = N * npara), diag(rep(TRUE, N * npara)))
    B_stack[non_zeros] <- B_par
  } else if (type == "countryidio") {
    non_zeros <- cbind(matrix(TRUE, ncol = p_joint, nrow = N * npara), ifelse(kronecker(diag(N), rep(1, npara)) == 1, TRUE, FALSE))
    B_stack[non_zeros] <- B_par
  } else if (type == "countryidio_nomu") {
    non_zeros <- cbind(matrix(TRUE, ncol = p_joint, nrow = N * npara), ifelse(kronecker(diag(N), c(rep(1, npara - 1), 0)) == 1, TRUE, FALSE))
    B_stack[non_zeros] <- B_par
  } else {
    stop("No valid model type selected.")
  }
  B_i <- makeBi(npara, N, p, p_joint, B_stack, type = type)
  
  return(list(B_stack = B_stack, B_i = B_i))
}
#' Loadings array
#'
#' @inheritParams Gibbs2_SM_SA_sampler
#' @inheritParams sumffkronV
#' @inheritParams bdiagByTime
#' @inheritParams makeBstart
#' @param B_stack describe later
#'
#' @return array of loadings
#' @export
makeBi <- function(npara, N, p, p_joint, B_stack, type = "allidio") { # type = c("allidio","countryidio")
  if (type == "allidio") {
    if (sum(B_stack) == 0) {
      B_i <- array(B_stack, dim = c(npara, npara + p_joint, N))
    } else {
      B <- aperm(array(c(t(B_stack)), dim = c(p_joint + npara * N, npara, N)), perm = c(2, 1, 3))
      B_i <- array(dim = c(npara, p_joint + npara, N))
      for (i in 1:N) {
        B_i[, , i] <- B[, colSums(B[, , i]) != 0, i]
      }
    }
  } else if (type == "countryidio" | type == "countryidio_nomu") {
    B <- aperm(array(c(t(B_stack)), dim = c(p_joint + N, npara, N)), perm = c(2, 1, 3))
    B_i <- array(dim = c(npara, p_joint + 1, N))
    for (i in 1:N) {
      B_i[, , i] <- B[, colSums(B[, , i]) != 0, i]
    }
  } else {
    stop("No valid model type selected.")
  }
  return(B_i = B_i)
}
#' Partial effects matrix and array
#'
#' @inheritParams sumffkronV
#' @inheritParams bdiagByTime
#' @param D_par number of D parameters
#'
#' @return partial effects matrix and array for starting
#' @export
makeDstart <- function(npara, N, nreg, D_par) {
  Dmat <- matrix(D_par, nrow = npara, ncol = nreg)
  Dstack <- kronecker(diag(N), Dmat)
  D_i <- array(D_par, c(npara, nreg, N))
  return(list(Dstack = Dstack, D_i = D_i))
}
set_scale_Vhat <- function(Vhat, inc_obs_old, inc_obs_new) {
  inc_obs_old / inc_obs_new * Vhat
}