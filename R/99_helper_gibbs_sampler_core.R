#' Sum in the inverse of the posterior variance for loadings/partial effects
#'
#' @param availableObs available observations
#' @param fPost backward sampled states (FFBS output)
#' @param w_regs regressor values, typically a matrix
#' @param Viarray VCOV array (npara x npara x TT) of cross-sectional unit i
#'
#' @return  cronecker summation part
#' @export
sumffkronV <- function(availableObs, fPost, w_regs, Viarray) {
  summ <- 0
  if (is.null(w_regs)) {
    for (tt in availableObs) {
      f <- fPost[, tt]
      V <- Viarray[, , tt]
      summ <- summ + kronecker(f %*% t(f), solve(V))
    }
  } else {
    for (tt in availableObs) {
      f <- c(fPost[, tt], w_regs[, tt])
      V <- Viarray[, , tt]
      summ <- summ + kronecker(f %*% t(f), solve(V))
    }
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
sumfyV <- function(availableObs, fPost, w_regs, yiObs, Viarray) {
  summ <- 0
  if (is.null(w_regs)) {
    for (tt in availableObs) {
      f <- fPost[, tt]
      y <- yiObs[, tt]
      V <- Viarray[, , tt]
      summ <- summ + solve(V) %*% y %*% t(f)
    }
  } else {
    for (tt in availableObs) {
      f <- c(fPost[, tt], w_regs[, tt])
      y <- yiObs[, tt]
      V <- Viarray[, , tt]
      summ <- summ + solve(V) %*% y %*% t(f)
    }
  }
  return(as.numeric(summ))
}
#' Matrix square root (by spectral decomposition)
#'
#' @param V Normal matrix
#'
#' @return matrix square root of V
#' @export
mat_sqrt <- function(V) {
  if (any(is.na(V))) {
    dim_V <- dim(V)[1]
    matrix(NA, nrow = dim_V, ncol = dim_V)
  } else {
    eig <- eigen(V)
    val <- eig$values
    vec <- eig$vectors
    # vec %*% diag(sqrt(val)) %*% t(vec)
    vec %*% diag(sqrt(val))
  }
}
#' Matrix square root of V-matrix
#'
#' Via helper [mat_sqrt()] using  spectral decomposition
#'
#' @param v_mat V-matrix
#' @param dim_v dimension of V-matrix
#' @param NN_times_TT number of cross sectional units (`NN`) multiplied by time
#'    periods (`TT`)
#'
#' @return matrix square root of V
#' @export
compute_mat_square_root_V <- function(v_mat, dim_v, NN_times_TT) {
  array(apply(v_mat, 3, mat_sqrt), c(dim_v, dim_v, NN_times_TT))
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
utuSum <- function(countryA, u, VhatSqrt, TT, NN, TT_num_y, npara) {
  uSplit <- lapply(
    split(u, matrix(rep(1:NN, each = TT_num_y), ncol = TT, byrow = T)),
    matrix, ncol = TT)
  sumUtu <- 0
  utildeSplit <- lapply(1:NN, matrix, data = NA, nrow = npara, ncol = TT)
  if (isTRUE(countryA)) {
    sumUtu_individ <- lapply(1:NN, matrix, data = NA, nrow = npara, ncol = npara)
  }
  for (i in 1:NN) {
    for (t in 1:TT) {
      utildeSplit[[i]][, t] <- solve(VhatSqrt[, , (i - 1) * TT + t]) %*% uSplit[[i]][, t]
    }
    utu <- apply(utildeSplit[[i]], 2, function(x) {
      x %*% t(x)
    })
    sumUtu_i <- apply(utu, 1, sum, na.rm = TRUE)
    sumUtu <- sumUtu + sumUtu_i

    if (isTRUE(countryA)) {
      sumUtu_i_mat <- matrix(sumUtu_i, ncol = npara)
      sumUtu_individ[[i]] <- (t(sumUtu_i_mat) + sumUtu_i_mat) / 2
    }
  }
  if (isTRUE(countryA)) return(sumUtu_individ)
  if (isFALSE(countryA)) return(matrix(sumUtu, ncol = npara))
  # countryA = TRUE -> $sumUtu_individ # countryA = FALSE -> $sumUtu_total
  # return(list(sumUtu_total = matrix(sumUtu, ncol = npara), sumUtu_individ = sumUtu_individ))
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
  # Dmat <- matrix(D_par, nrow = npara, ncol = nreg)
  # Dstack <- kronecker(diag(N), Dmat)
  # D_i <- array(D_par, c(npara, nreg, N))
  # return(list(Dstack = Dstack, D_i = D_i))
  D_i <- array(D_par, c(npara, nreg, N))
  Dstack <- as.matrix(Matrix::bdiag(plyr::alply(D_i, 3)))

  return(list(Dstack = Dstack, D_i = D_i))
}
set_scale_Vhat <- function(Vhat, inc_obs_old, inc_obs_new) {
  inc_obs_old / inc_obs_new * Vhat
}
set_f_out <- function(num_fac, TT, itermax) {
  array(0, dim = c(num_fac, TT, itermax))
}
# speichert den Output von FFBS
# cSTORE <- matrix(rep(0,N_num_y * itermax), ncol = itermax)
set_B_out <- function(N_num_y, num_fac, itermax) {
  array(0, dim = c(N_num_y, num_fac, itermax))
}
# speichert die Gibbs-Zuege der Ladungen auf die latenten Faktoren
# uSTORE <- array(0, dim = c(N_num_y, TT, itermax) )
set_D_out <- function(N_num_y, NN, nreg, itermax) {
  if (nreg != 0) {
    # speichert die Koeffizienten der Makroregressor-Koeffizienten
    return(list(DSTORE = array(0, dim = c(N_num_y, NN * nreg, itermax)),
                DiSTORE = array(0, dim = c(N_num_y / NN, nreg, NN))))
  } else {
    return(NULL)
  }
}
set_A_out <- function(countryA, num_y, itermax, NN){
  if (countryA) {
    array(0, dim = c(num_y, num_y, itermax, NN))
  } else {
    array(0, dim = c(num_y, num_y, itermax))
  }
}
set_V_out <- function(sampleV, n_num_y, num_mcmc) {
  if (isTRUE(sampleV)) return(NULL)
  matrix(0, nrow = n_num_y, ncol = num_mcmc)
}
set_V_tmp <- function(V_DIAG_EST, V_HAT_DIAG_SCALE,
                      Vhat, VhatDiagScale_start,
                      NN, TT, NN_TT, N_num_y, num_y) {
  out_list <- list(
    Vstart = NULL,
    Vhat = Vhat,
    VhatFix = NULL,
    VhatArrayBdiagByTimeFix = NULL
  )
  if (isFALSE(V_DIAG_EST)) return(out_list)
  if (isTRUE(V_DIAG_EST)) {
    # Vstart <- initials$Vhat[1, 1, 1]
    out_list$Vstart <- Vhat[1, 1, 1]
    if (V_HAT_DIAG_SCALE) {
      # speichert Vhat Initialisierung
      out_list$VhatFix <- Vhat
      # Sortiert Vhat nach der Zeit: wie Vhat eine num_y x num_y x N*TT Array,
      # aber die Elemente [1:num_y, 1:num_y, 1:N] sind nun die
      # Messfehlervarianzen aller Laender zum ersten Zeitpunkt,
      #[1:num_y, 1:num_y, (N+1):(2N)] alle Laender zum zweiten Zeitpunkt, usw.
      out_list$VhatArrayBdiagByTimeFix <- bdiagByTime(out_list$VhatFix,
                                                      npara = num_y,
                                                      N = NN, TT = TT,
                                                      Nnpara = N_num_y)
      if (!missing(VhatDiagScale_start) && !is.null(VhatDiagScale_start)) {
        out_list$Vhat <- array(
          sapply(1:(NN_TT), \(x) VhatDiagScale_start[, , x] %*% Vhat[, , x]),
          dim = c(num_y, num_y, NN_TT)
        ) # setzt Vhat im Fall VhatDiagScale = TRUE
        out_list$Vstart <- VhatDiagScale_start[1, 1, 1]
      } else {
        out_list$Vstart <- 1
      }
    }
  }
  return(out_list)
}
set_A_country_array <- function(countryA, A, NN, num_y) {
  if (isTRUE(countryA)) {
    A_country <- plyr::alply(replicate(NN, A), 3) # Liste: enthaelt fuer jedes Land eine eigene A Matrix nach dem Vorbild des uebergebenen A's.
    array(unlist(A_country), c(num_y, num_y, NN)) # macht aus der Liste ein Array
  } else if (isFALSE(countryA)) {
    NULL
  }
}
compute_V_hat_array_A <- function(VhatSqrt, A,
                                  countryA, A_countryArray,
                                  NN, TT, NN_TT, num_y) {
  if (countryA) {
    VhatArray_A <- array(0, dim = c(num_y, num_y, NN_TT))
    for (i in 1:NN) {
      VhatArray_A[, , (TT * (i - 1) + 1):(i * TT)] <- array(
        apply(VhatSqrt[, , (TT * (i - 1) + 1):(i * TT)],
              3, \(x) x %*% A_countryArray[, , i] %*% t(x)),
        c(num_y, num_y, TT)
      )
    }
    return(VhatArray_A)
  } else {
    array(
      apply(VhatSqrt, 3, function(X) {X %*% A %*% t(X)}),
      c(num_y, num_y, NN_TT)
    )
  }
}
sample_V <- function(VhatDiagScale, VhatArrayBdiagByTimeFix, VhatFix, u,
                     TT, NN_TT, num_y, availableObs_crossSection,
                     alpha0, beta0) {
  if (VhatDiagScale) {
    u <- sapply(1:TT, \(x) (1 / sqrt(diag(VhatArrayBdiagByTimeFix[, , x]))) * u[, x])
  }
  ssu <- apply(u, 1, \(x) sum(x^2, na.rm = TRUE))
  TT_available <- rep(apply(availableObs_crossSection, 1, sum), each = num_y)
  ssu_TT <- cbind(ssu, TT_available)
  Vdiag <- apply(ssu_TT, 1, \(x) LaplacesDemon::rinvgamma(1, shape = alpha0 + x[2] / 2, scale = beta0 + x[1] / 2))
  VdiagMat <- matrix(Vdiag, ncol = num_y, byrow = T)
  VdiagMat_extend <- VdiagMat[rep(1:nrow(VdiagMat), each = TT), ]
  VdiagArray <- array(apply(VdiagMat_extend, 1, diag), dim = c(num_y, num_y, NN_TT))
  if (VhatDiagScale) {
    VdiagArray <- array(sapply(1:(NN_TT), \(x) VdiagArray[, , x] %*% VhatFix[, , x]), dim = c(num_y, num_y, NN_TT))
  }
  VhatSqrt <- compute_mat_square_root_V(VdiagArray, num_y, NN_TT)
  return(list(V_diag = Vdiag, V_hat_sqrt = VhatSqrt))
}
compute_residuals <- function(y, f, B, D, regs) {
  y - B %*% f - D %*% regs
}
sample_A <- function(countryA, diagA, scaleA,
                     u, VhatSqrt,
                     NN, TT, num_y, TT_num_y,
                     NN_TT_avail, availableObs_crossSection,
                     shape0, rate0, Psi0, nu0) {
  if (countryA) {
    utu_country <- utuSum(countryA = countryA, u = u, VhatSqrt = VhatSqrt, NN = NN, TT = TT, TT_num_y = TT_num_y, npara = num_y)
    if (diagA) {
      uSum_country <- lapply(utu_country, diag)
      Adiag_country <- lapply(1:NN, \(xx)  sapply(uSum_country[[xx]], \(x) LaplacesDemon::rinvgamma(n = 1, shape = shape0 + 0.5 * sum(availableObs_crossSection[xx, ]), scale = rate0 + 0.5 * x)))
      A_countryArray <- array(unlist(lapply(Adiag_country, diag)), c(num_y, num_y, NN))
    } else {
      A_country <- lapply(1:NN, \(x) LaplacesDemon::rinvwishart(nu = sum(availableObs_crossSection[x, ]) + nu0, S = Psi0 + utu_country[[x]]))
      A_countryArray <- array(unlist(A_country), c(num_y, num_y, NN))
    }
    return(A_countryArray)
  } else {
    utu <- utuSum(countryA = countryA, u = u, VhatSqrt = VhatSqrt, NN = NN, TT = TT, TT_num_y = TT_num_y, npara = num_y)
    utu <- (utu + t(utu)) / 2
    if (diagA) {
      # uSplitSqrt <- lapply(uSplit,\(x) x^2)
      # uSum <- apply(sapply(uSplitSqrt, \(x) apply(x,1,sum)),1,sum)
      uSum <- diag(utu)
      A_diag <- sapply(uSum, \(x) invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (NN_TT_avail), rate = rate0 + 0.5 * x))
      A <- diag(A_diag)
    } else if (scaleA) {
      uSum <- sum(diag(utu))
      A_scale <- invgamma::rinvgamma(n = 1, shape = shape0 + 0.5 * (num_y * NN_TT_avail), rate = rate0 + 0.5 * uSum)
      A <- diag(rep(A_scale, num_y))
    } else {
      A <- LaplacesDemon::rinvwishart(nu = NN_TT_avail + nu0, S = Psi0 + utu)
    }
    return(A)
  }
}
compute_B_post_full <- function(availableObs,
                                invOmega0, B0_D0,
                                fPost, w_regs,
                                Viarray, yiObs, selectR,
                                try_catch_errors) {
  B_Omega <- compute_Omega1(invOmega0,
                            availableObs,
                            selectR,
                            fPost,
                            w_regs,
                            Viarray,
                            try_catch_errors)
  B_mean <- compute_B_mean(B_Omega,
                           invOmega0,
                           B0_D0,
                           availableObs, selectR,
                           fPost, w_regs,
                           yiObs,  Viarray)
  B_Sigma <- compute_Sigma_adjust(B_Omega)
  return(list(B_mean = B_mean, B_Sigma = B_Sigma))
}
compute_B_mean <- function(Omega,
                           invOmega0,
                           B0_D0,
                           availableObs, selectR,
                           fPost, w_regs,
                           yiObs,  Viarray) {
  beta1_mid <- sumfyV(availableObs,
                      fPost = fPost,
                      w_regs = w_regs,
                      yiObs = yiObs,
                      Viarray = Viarray)
  # drop(Omega %*% (selectR %*% (beta1_mid + invOmega_B0_D0)))
  drop(Omega %*% (selectR %*% beta1_mid + invOmega0 %*% selectR %*% B0_D0))
  # return(Omega %*% (selectR %*% (c(beta1_mid) + invOmega %*% c(B0, D0))))
  # beta1 <- Omega1 %*% (c(beta1_mid) + invOmega0 %*%  c(B0[,,i]) )
  # beta1 <- Omega1 %*% (selectR %*% (c(beta1_mid) + invOmega0 %*% c(B0[, , i], D0[, , i])))
}
compute_Omega1 <- function(invOmega0 ,
                           availableObs,
                           selectR,
                           fPost,
                           w_regs,
                           Viarray,
                           try_catch_errors) {
  ## Posterior Momente fuer die Ladungen und die partiellen Effekte
  invOmega1_part2 <- sumffkronV(availableObs,
                                fPost = fPost,
                                w_regs = w_regs,
                                Viarray = Viarray)

  invOmega1 <- invOmega0 + selectR %*% invOmega1_part2  %*% t(selectR)
  # invOmega1 <- invOmega0 + invOmega1_part2

  if (isFALSE(is.null(try_catch_errors))) {
    tryCatch(
      {
        solve(invOmega1)
        # solve(selectR %*% invOmega1 %*% t(selectR))
      },
      error = function(e) {
        # print(invOmega1)
        if (try_catch_errors$storePath != "none") {
          if (try_catch_errors$store_count == 0) {
            dir.create(try_catch_errors$storePath_adj, recursive = TRUE)
            # store_count <- try_catch_errors$store_count + 1
          }
          saveRDS(
            invOmega1,
            file = paste0(try_catch_errors$storePath_omg, "_",
                          try_catch_errors$iter))
        }
        solve(invOmega1, tol = 0)
        # solve(selectR %*% invOmega1 %*% t(selectR), tol = 0)
      }
    )
  }
  return(solve(invOmega1))
  # return(solve(selectR %*% invOmega1 %*% t(selectR)))
}
compute_Sigma_adjust <- function(Omega) {
  0.5 * Omega + 0.5 * t(Omega)
  # Sigma <- 0.5 * (selectR %*% Omega1 %*% t(selectR)) + 0.5 * t(selectR %*% Omega1 %*% t(selectR))
}
get_identificiation_restrictions <- function(type, num_joint_fac,
                                             num_y, num_reg) {
  ## Identifikationsrestriktionen fuer die Ladungen
  # Dselect <- diag(6)
  if (type == "allidio") {
    if (num_joint_fac == 0) {
      # muss fuer gemeinsamen faktor angepasst werden
      lower <- c(rep(0, num_y), rep(-Inf, num_y * num_reg))
      # muss fuer gemeinsamen faktor angepasst werden
      upper <- rep(Inf, num_y + num_y * num_reg)
    } else if (num_joint_fac == 1) {
      # lower bound for factors
      lower <- c(rep(0, num_y), rep(-Inf, num_y * num_reg))
      lower_first <- c(0, rep(-Inf, num_y - 1), lower)
      lower_rest  <- c(rep(-Inf, num_y), lower)
      lower <- rbind(lower_first, lower_rest)
      # upper bound for factors
      upper <- rep(Inf, 2 * num_y + num_y * num_reg)
      upper <- rbind(upper, upper)
      # lower_first <- c(0, rep() rep(-Inf, 2), rep(0, 3))
      # if (i == 1) {
      #   #   Dselect <- Dselect[-c(2,3),]
      #   lower <- c(0, rep(-Inf, 2), rep(0, 3))
      # } else {
      #   # Dselect <- Dselect[-c(1,2,3),]
      #   lower <- c(rep(-Inf, 3), rep(0, 3))
      # }
      # upper <- rep(Inf, 6)
    } else {
      stop("Not yet implemented.")
    }
  } else if (type == "countryidio") {
    if (i == 1) {
      lower <- c(0, rep(-Inf, 2), 0, rep(-Inf, 2))
      #   Dselect <- Dselect[-c(2,3),]
    } else {
      lower <- c(rep(-Inf, 3), 0, rep(-Inf, 2))
      #   Dselect <- Dselect[-c(1,2,3),]
    }
    upper <- rep(Inf, 6)
  } else if (type == "countryidio_nomu") {
    if (i == 1) {
      lower <- c(0, rep(-Inf, 2), 0, -Inf)
      #   Dselect <- Dselect[-c(2,3),]
    } else {
      lower <- c(rep(-Inf, 3), 0, -Inf)
      #   Dselect <- Dselect[-c(1,2,3),]
    }
    upper <- rep(Inf, 5)
  } else {
    stop("No valid model type selected.")
  }
  return(list(upper = upper, lower = lower))
}
sample_B_D <- function(mean_B_full, sigma_B_full, upper, lower, num_jnt_fac, num_y, LEGACY = FALSE) {
  # BDsamp <- tmvnsim::tmvnsim(1, length(upper), lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma)$samp
  # Bvec <- MASS::mvrnorm(n = 1, mu = selectR %*% beta1 + selectC, Sigma = selectR %*% Omega1 %*% t(selectR))
  # Bvec <- as.numeric(tmvtnorm::rtmvnorm(n = 1, mean = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma,
  #                          lower = lower, upper = upper))

  # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(selectR %*% beta1 + selectC), sigma = Sigma )$samp
  # Bvec <- tmvnsim::tmvnsim(1,length(upper),lower = lower, upper = upper, means = as.numeric(beta1 + selectC), sigma = Sigma )$samp
  # browser()
  if (isTRUE(LEGACY)) {
    Bsamp_full <- tmvnsim::tmvnsim(1, length(upper), lower = lower, upper = upper,
                                   means = mean_B_full, sigma = sigma_B_full)$samp

  } else {
    # Bsamp_full <- tmvtnorm::rtmvnorm(1, mean = mean_B_full,
    #                                  sigma = sigma_B_full,
    #                                  lower = lower, upper = upper)
    Bsamp_full <- tmvtnsim:::rtmvnormcpp(mean = matrix(mean_B_full, ncol = ncol(sigma_B_full)),
                                         sigma = sigma_B_full,
                                         blc = diag(ncol(sigma_B_full)),
                                         lower = matrix(lower, ncol = ncol(sigma_B_full)),
                                         upper = matrix(upper, ncol = ncol(sigma_B_full)),
                                                        init = matrix(0, ncol = ncol(sigma_B_full)),
                                                        burn = 10)
  }

  get_BD_samp(Bsamp_full, num_jnt_fac, num_y)
}
get_BD_samp <- function(Bsamp_full, num_jnt_fac, num_y) {
  Bvec <- Bsamp_full[1:((num_jnt_fac + 1) * num_y)]
  if (num_jnt_fac != 0) {
    idioFac <- Bvec[-(1:(num_y * num_jnt_fac))]
    if (type == "countryidio_nomu") {
      idioFac <- c(idioFac, 0)
    }
    # idioPos <- TRUE
    jointFac <- matrix(Bvec[1:(num_y * num_jnt_fac)], ncol = num_jnt_fac)
    # neu
    # Unterscheidung i =1 und sonst
    # jointFac <- matrix(0, nrow = num_y, ncol = num_jnt_fac)
    # jointFac_select <- lower.tri(jointFac, diag = T)
    # jointFac[jointFac_select] <- Bvec[1:sum(jointFac_select)]
    # idioFac <- Bvec[-(1:sum(jointFac_select))]
  } else {
    jointFac <- NULL
    idioFac <- Bvec # muss fuer regs geaendert werden bzw. Bvec oben neu definieren (D abziehen)
  }
  # idioPos <- all(idioFac > 0)
  # gewekePos <- TRUE

  # if(i == 1 & num_jnt_fac != 0){
  #   gewekeMatrix <- jointFac[1:num_jnt_fac, 1:num_jnt_fac] # muss verallgemeinert werden fuer num_y < num_jnt_fac
  #   gewekePos <- all(diag(gewekeMatrix, num_jnt_fac) > 0) # muss verallgemeinert werden fuer num_y < num_jnt_fac
  #   #gewekePos <- TRUE
  # }
  # if(gewekePos & idioPos){valid <- TRUE}
  # if(!identification){valid <- TRUE}
  # valid <- TRUE
  # ident_control <- ident_control + 1
  # }

  # if(ident_block){break}
  return(list(Bsamp_idi = idioFac, Bsamp_jnt = jointFac,
              Dsamp = Bsamp_full[-(1:((num_jnt_fac + 1) * num_y))]))
}
compute_invOmega0_B0_D0 <- function(invOmega0, B0, D0) {
  NN <- dim(B0)[3]
  out <- list(NN)
  for (n in seq_len(NN)) {
    B0_D0 <- c(B0[, , n], D0[, , n])
    out[[n]] <- drop(invOmega0 %*% B0_D0)
  }
  return(out)
}
get_B <- function(Bjoint, Bidio, num_fac_jnt, NN, num_y, type) {
  if (type == "allidio") {
    if (num_fac_jnt != 0) {
      B <- cbind(Bjoint, diag(c(Bidio)))
    } else {
      B <- diag(c(Bidio))
    }
  } else if (type == "countryidio" | type == "countryidio_nomu") {
    BidioMat <- as.matrix(Matrix::bdiag(plyr::alply(array(Bidio, c(num_y, 1, NN)), 3)))
    if (num_fac_jnt != 0) {
      B <- cbind(Bjoint, BidioMat)
    } else {
      B <- BidioMat
    }
  }
  return(B)
}
get_id_fpost <- function(num_fac_jnt, num_y, NN, type) {
  id_out <- matrix(0, ncol = NN, nrow = num_fac_jnt + num_y)
  for (i in 1:NN) {
    if (type == "allidio") {
      if (num_fac_jnt != 0) {
        # generalize of no num_fac_jnt is included
        id_tmp <- c(1:num_fac_jnt, (num_fac_jnt + 1 + num_y * (i - 1)):(num_fac_jnt + num_y * i))
      } else {
        id_tmp <- (num_fac_jnt + 1 + num_y * (i - 1)):(num_fac_jnt + num_y * i)
      }
    } else if (type == "countryidio" | type == "countryidio_nomu") {
      if (num_fac_jnt != 0) {
        id_tmp <- c(1:num_fac_jnt, num_fac_jnt + i) # generalize of no num_fac_jnt is included
      } else {
        id_tmp <- i
      }
    } else {
      stop("No valid model type selected.")
    }
    id_out[, i] <- id_tmp
  }
  return(id_out)
}
get_id_wreg <- function(num_w_reg, NN) {
  out_id <- matrix(0, nrow = num_w_reg, NN)
  for (i in seq_len(NN)) {
    out_id[, i] <- (1 + num_w_reg * (i - 1)):(i * num_w_reg)
  }
  out_id
}
sample_B_full_cpp <- function(yObs, availableObs_crossSection,
                              fPost, VhatArray_A, w_reg_info,
                              Omega0, B0, D0,
                              id_f, selectR, lower, upper,
                              NN, TT, N_num_y, num_y, num_fac_jnt, num_par_all,
                              type, B_trues, D_trues) {
  # Ergebnis-Matrix fuer die gemeinsamen Ladungen (wird spaeter befuellt)
  Bjoint <- matrix(rep(0, N_num_y * num_fac_jnt), ncol = num_fac_jnt)
  # Ergebnis-Matrix fuer die idiosynkratischen Ladungen (wird spaeter befuellt)
  Bidio   <- matrix(rep(0, N_num_y), ncol = NN)
  nregs   <- w_reg_info$nregs
  if(nregs == 0){
    Dregs <- NULL
    Dregs_i <- NULL
  } else {
    Dregs   <- matrix(0, N_num_y, NN * nregs)
    Dregs_i <- array(0, dim = c(num_y, nregs, NN))
  }
  invOmega0_i <- solve(Omega0)
  upper_taken <- upper
  lower_taken <- lower
  for (i in 1:NN) {
    if (num_fac_jnt > 0) {
      bound_id <- 2
      if (i == 1) {
        bound_id <- 1
      }
      upper_taken <- upper[bound_id, ]
      lower_taken <- lower[bound_id, ]
    }
    # Vhat Array for cross-sectional unit (country) i
    Viarray <- VhatArray_A[, , (1 + (i - 1) * TT):(i * TT)]
    # yObs for cross-sectional unit (country) i
    yiObs <- yObs[(1 + num_y * (i - 1)):(num_y * i), ]
    # availableObs <- which(!is.na(yObs[1 + num_y * (i-1),]))
    availableObs <- which(availableObs_crossSection[i, ])
    f_post_i <- fPost[id_f[, i], ]
    # invOmega0_B0_D0_i <- invOmega0_B0_D0[[i]]
    B0_D0_i  <- c(B0[, , i], D0[, , i])
    if(nregs == 0) {
      w_regs_i <- matrix(0, 0, 0)
    } else {
      w_regs_i <- w_reg_info$w_reg[w_reg_info$id_reg[, i], ]
    }

    B_D_samp <-  sample_B_D_cpp(availableObs, invOmega0_i, B0_D0_i,
                                f_post_i, w_regs_i, Viarray, yiObs,
                                selectR,
                                upper = upper_taken,
                                lower = lower_taken,
                                num_fac_jnt, num_y)

    Dsamp <- B_D_samp$Dsamp
    Bsamp_jnt <- B_D_samp$Bsamp_jnt
    Bsamp_idi <- B_D_samp$Bsamp_idi
    # Bsamp_jnt <- B_trues[(i - 1) * 3 + 1:3, 1]
    # Bsamp_idi <- diag(B_trues[, 1:30])[(i - 1) * 3 + 1:3]
    # Bsamp_idi <- diag(B_trues[, 2:31])[(i - 1) * 3 + 1:3]
    if (nregs != 0) {
      # Dsamp <- D_trues[(1 + num_y * (i - 1)):(i * num_y),
      #       (1 + nregs * (i - 1)):(i * nregs)]
      Dregs[(1 + num_y * (i - 1)):(i * num_y),
            (1 + nregs * (i - 1)):(i * nregs)] <- Dsamp
      Dregs_i[, , i] <- Dsamp
    }
    if (num_fac_jnt != 0) {
      Bjoint[((i - 1) * num_y + 1):(i * num_y), ] <- Bsamp_jnt
    }
    Bidio[, i] <- Bsamp_idi
  }
  Bfacs <- get_B(Bjoint, Bidio, num_fac_jnt = num_fac_jnt,
                 NN = NN, num_y = num_y, type = type)
  Bfacs_i <- makeBi(npara = num_y, N = NN, p = num_par_all,
                    p_joint = num_fac_jnt, B_stack = Bfacs,
                    type = type)
  return(list(Bfacs = Bfacs, Dregs = Dregs,
              Bfacs_i = Bfacs_i, Dregs_i = Dregs_i))

}
sample_B_full <- function(yObs, availableObs_crossSection,
                          fPost, VhatArray_A, w_reg_info,
                          Omega0, B0, D0,
                          id_f, selectR, lower, upper,
                          NN, TT, N_num_y, num_y, num_fac_jnt, num_par_all,
                          type) {
  # Ergebnis-Matrix fuer die gemeinsamen Ladungen (wird spaeter befuellt)
  Bjoint <- matrix(rep(0, N_num_y * num_fac_jnt), ncol = num_fac_jnt)
  # Ergebnis-Matrix fuer die idiosynkratischen Ladungen (wird spaeter befuellt)
  Bidio   <- matrix(rep(0, N_num_y), ncol = NN)
  nregs   <- w_reg_info$nregs
  Dregs   <- matrix(0, N_num_y, NN * nregs)
  Dregs_i <- array(0, dim = c(num_y, nregs, NN))
  invOmega0_i <- solve(Omega0)
  upper_taken <- upper
  lower_taken <- lower
  for (i in 1:NN) {
    if (num_fac_jnt > 0) {
      bound_id <- min(i, 2)
      upper_taken <- upper[bound_id, ]
      lower_taken <- lower[bound_id, ]
    }
    # Vhat Array for cross-sectional unit (country) i
    Viarray <- VhatArray_A[, , (1 + (i - 1) * TT):(i * TT)]
    # yObs for cross-sectional unit (country) i
    yiObs <- yObs[(1 + num_y * (i - 1)):(num_y * i), ]
    # availableObs <- which(!is.na(yObs[1 + num_y * (i-1),]))
    availableObs <- which(availableObs_crossSection[i, ])
    f_post_i <- fPost[id_f[, i], ]
    # invOmega0_B0_D0_i <- invOmega0_B0_D0[[i]]
    B0_D0_i  <- c(B0[, , i], D0[, , i])
    if(nregs == 0){
      w_regs_i <- NULL
    } else {
      w_regs_i <- w_reg_info$w_reg[w_reg_info$id_reg[, i], ]
    }

    B_post <- compute_B_post_full(availableObs = availableObs,
                                  invOmega0 = invOmega0_i,
                                  B0_D0 = B0_D0_i,
                                  fPost = f_post_i,
                                  w_regs = w_regs_i,
                                  Viarray = Viarray,
                                  yiObs = yiObs,
                                  selectR = selectR,
                                  try_catch_errors = NULL)
    # try_catch_errors = list(storePath = storePath,
    #                         storePath_adj = storePath_adj,
    #                         storePath_omg = storePath_omg,
    #                         store_count = store_count,
    #                         iter = iter)
    # )
    bmean <- B_post$B_mean
    Sigma <- B_post$B_Sigma
    B_D_samp <- sample_B_D(mean_B_full = bmean, Sigma,
                           upper = upper_taken,
                           lower = lower_taken,
                           num_jnt_fac = num_fac_jnt,
                           num_y = num_y)
    Dsamp <- B_D_samp$Dsamp
    Bsamp_jnt <- B_D_samp$Bsamp_jnt
    Bsamp_idi <- B_D_samp$Bsamp_idi
    if (nregs != 0) {
      Dregs[(1 + num_y * (i - 1)):(i * num_y),
            (1 + nregs * (i - 1)):(i * nregs)] <- Dsamp
      Dregs_i[, , i] <- Dsamp
    }
    if (num_fac_jnt != 0) {
      Bjoint[((i - 1) * num_y + 1):(i * num_y), ] <- Bsamp_jnt
    }
    Bidio[, i] <- Bsamp_idi
  }
  Bfacs <- get_B(Bjoint, Bidio, num_fac_jnt = num_fac_jnt,
                 NN = NN, num_y = num_y, type = type)
  Bfacs_i <- makeBi(npara = num_y, N = NN, p = num_par_all,
                    p_joint = num_fac_jnt, B_stack = Bfacs,
                    type = type)
  return(list(Bfacs = Bfacs, Dregs = Dregs,
              Bfacs_i = Bfacs_i, Dregs_i = Dregs_i))

}
get_outlist_Gibbs_sampler <- function(fSTORE, BSTORE, DSTORE,
                                      ASTORE, VSTORE, uSTORE,
                                      BD0STORE, Omega0STORE,
                                      block_count,
                                      msg_error_kf,
                                      initials) {
  out_list <- list(f = fSTORE, B = BSTORE, D = DSTORE,
                   A = ASTORE, V = VSTORE, u = uSTORE,
                   BD0STORE = BD0STORE,
                   Omega0STORE = Omega0STORE,
                   blockCount = block_count,
                   errorMsg = msg_error_kf,
                   initials = initials)
  return(out_list)
}