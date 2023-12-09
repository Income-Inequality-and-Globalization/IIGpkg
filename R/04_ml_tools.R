#' Helper to compute V-hat adjusted matrix to the dimension of Y.
#'
#' @param V_sq_adj 3x3 or other dimensional V-hat matrix
#' @param NN_TT number of cross sectional units `x` time periods
#' @param num_y dimension of the output V-hat matrix (should be smaller then the
#'    current dimension of the `V_sq_adj` as passed via the first argument)
#'
#' @return adjusted matrix
#' @export
get_V_sq_adj2 <- function(V_sq_adj, NN_TT, num_y, num_v_sq_taken) {
  if (num_y > 1) {
    V_sq_adj2  <- replicate(NN_TT, V_sq_adj[1:num_y, 1:num_y, num_v_sq_taken]) # V_sq_adj[, , 92] or V_sq_adj[, , 199]
  } else {
    V_sq_adj2  <- array(replicate(NN_TT, V_sq_adj[1:num_y, 1:num_y, num_v_sq_taken]),
                        dim = c(1, 1, NN_TT))# V_sq_adj[, , 92] or V_sq_adj[, , 199]
  }
  return(V_sq_adj2)
}
#' Generate errors for Y process simulation
#'
#' @inheritParams get_V_sq_adj2
#' @param A_means 
#' @param num_y 
#' @param TT number of time periods
#' @param NN number of cross sectional units
#'
#' @return simulated errors
#' @export
get_u_simul <- function(V_sq_adj, A_means, num_y, TT, NN) {
  it <- 1
  for (ii in 1:NN) {
    for (tt in 1:TT) {
      Sigma_tmp <- as.matrix(V_sq_adj2[, , ((ii - 1) * TT + tt)]) %*% as.matrix(A_means[1:num_y, 1:num_y]) %*% t(as.matrix(V_sq_adj2[, , ((ii - 1) * TT + tt)]))
      # Sigma_tmp <- V_sq_adj[, , 120] %*% A_means %*% t(V_sq_adj[, , 120])
      Sigma_true[ , , it] <- Sigma_tmp
      if (is.na(Sigma_tmp[1, 1])) {
        u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- rep(NA, times = num_y)
      } else {
        u_simul[((ii - 1) * num_y + 1):(ii * num_y), tt] <- MASS::mvrnorm(
          1,
          rep(0, times = num_y),
          Sigma_true[ , , it])
      }
      it <- it + 1
    }
  }
  return(u_simul)
}
#' Simulate random walk factors
#'
#' @inheritParams get_V_sq_adj2
#' @inheritParams get_u_simul
#' @param num_joint_facs number of joint factors
#'
#' @return simulated random walk factors
#' @export
generate_f_random_walk <- function(NN_TT, num_y, TT, num_joint_facs) {
  f_simul <- matrix(rnorm(NN_TT * num_y + TT * njfacs, mean = 0, sd = 1),
                    nrow = NN * num_y + njfacs, ncol = TT)
  f_simul <- t(apply(f_simul, 1, cumsum)) 
  return(f_simul)
}
#' Simulate measurements/data/observations
#'
#' @param u errors
#' @param B loadings
#' @param f factors
#' @param D marginal regressor effects
#' @param wRegs regressors
#'
#' @return simulated data/measurements/observations
#' @export
generate_y_simul <- function(u, B, f, D = NULL, wRegs = NULL) {
  if (is.null(D) && is.null(wRegs)) return(B %*% f + u)
  if (is.null(D) || is.null(wRegs)) stop("Both, D and wRegs must NULL.")
  B %*% f + D %*% wRegs + u
}
#' Helper function to plat
#'
#' @param yObs simulated data/measurements/observations as produced via
#'    [IIGpkg::generate_y_simul]
#' @param f_simul simulated factors as produced via 
#'    [IIGpkg::generate_f_random_walk]
#' @inheritParams generate_f_random_walk
#' @inheritParams get_u_simul
#'
#' @return pure side effect function for plotting
#' @export
plot_data_ml_simul <- function(yObs, f_simul, num_y, NN, num_joint_facs) {
  y_max <- max(yObs)
  y_min <- min(yObs)
  par(mfrow = c(2, 1))
  plot(yObs[1, ], type = "l", ylim = c(y_min, y_max),
       ylab = "measurements", xlab = "time", main = "Measurements")
  for (tt in 2:(num_y * NN)) {
    lines(yObs[tt, ], type = "l")
  }
  y_max <- max(f_simul)
  y_min <- min(f_simul)
  plot(f_simul[1, ], type = "l", ylim = c(y_min, y_max), col = "blue",
       ylab = "factors", xlab = "time", main = "Factors -- blue: common factor")
  for (tt in 2:(num_y * NN + num_joint_facs)) {
    lines(f_simul[tt, ], type = "l")
  }
  par(mfrow = c(1, 1))
  return(invisible(NULL))
}
#' ML objective function for first identification scheme
#' 
#' This identification scheme fixes the KF-transition variances to 1.
#'
#' @param x argument for parameter vector over which to maximize
#' @param yObs observations/measurements/data
#' @param wReg regressors
#' @param B_means_used true values of loadings
#' @param A_means_used true values of adjustment matrix
#' @param Q_used true values of transition variances
#' @param V_sq_adj VCM matrix of measurements/data/observations
#' @param Phi random walk matrix
#' @param ccc replacmenet for (non-existent) constant
#' @param initX initialization value of X
#' @param initU initialization value of U
#' @param initP initialization value of P 
#' @param id_A_mat_optim integer vector giving the id's of the adjustment matrix
#'    over which to maximize
#' @inheritParams get_V_sq_adj2
#' @inheritParams get_u_simul
#'
#' @return log-likelihood value
#' @export
ml_objective_identification_01 <- function(x, yObs, wReg,
                                           B_means_used,
                                           A_means_used,
                                           Q_used, V_sq_adj,
                                           Phi, ccc,
                                           initX, initU, initP,
                                           id_A_mat_optim,
                                           NN, TT, num_y, NN_TT,
                                           num_jnt_fac, PLOT = FALSE) {
  B_means_used[, 1] <- x[1:(NN * num_y)]
  diag(B_means_used[, (num_jnt_fac + 1):ncol(B_means_used)]) <- x[(NN * num_y + 1):(2 * NN * num_y)]
  A_means_used[id_A_mat_optim] <- x[(2 * NN_num_y + 1):(2 * NN_num_y + num_a_optim)]
  VhatArray_A <- array(
        apply(V_sq_adj[1:num_y, 1:num_y, , drop = FALSE], 3,
              function(X) {
                X %*% as.matrix(A_means_used[1:num_y, 1:num_y]) %*% t(X)
                }
        ),
        c(num_y, num_y, NN_TT)
  )
  VhatArrayBdiagByTime <- IIGpkg::bdiagByTime(VhatArray_A, num_y, NN, TT, NN_num_y)

  out <- RcppSMCkalman::kfLGSSM(yObs = yObs,
                                uReg = NULL,
                                wReg = wReg,
                                A = Phi,
                                B = NULL,
                                C = B_means_used,
                                D = ccc,
                                Q = Q_used,
                                R = VhatArrayBdiagByTime,
                                initX = initX,
                                initU = initU,
                                initP = initP,
                                computeMFD = FALSE,
                                computeMSD = FALSE,
                                computePDD = TRUE,
                                computeLLH = TRUE)$kfLLHout
  if (out == Inf || is.na(out)) out <- -Inf
  if (PLOT) return(out)
  return(-out)
}
#' ML objective function for second identification scheme
#' 
#' This identification scheme fixes the idiosyncratic loadings to 1 and first 
#' joint loading to 1.
#'
#' @inheritParams ml_objective_identification_01
#'
#' @return log-likelihood value
#' @export
ml_objective_identification_02 <- function(x, yObs, wReg,
                                           B_means_used,
                                           A_means_used,
                                           Q_used, V_sq_adj,
                                           Phi, ccc,
                                           initX, initU, initP,
                                           id_A_mat_optim,
                                           NN, TT, num_y, NN_TT,
                                           num_jnt_fac, PLOT = FALSE) {
  id_B_optim_fr <- seq(
    from = 1,
    to = (NN_num_y * num_jnt_fac - 1),
    by = 1)
  id_Q_optim_fr <- seq(
    from = (NN_num_y * num_jnt_fac),
    to = (NN_num_y * num_jnt_fac + NN_num_y + num_jnt_fac - 1),
    by = 1)
  id_A_optim_fr <- seq(
    from = (NN_num_y * num_jnt_fac + NN_num_y + num_jnt_fac),
    to = (NN_num_y * num_jnt_fac + NN_num_y + num_jnt_fac + length(id_A_mat_optim) - 1),
    by = 1)

  B_means_used[2:(NN_num_y), 1] <- x[id_B_optim_fr]
  diag(Q_used)                  <- (x[id_Q_optim_fr])^2
  A_means_used[id_A_mat_optim]  <- (x[id_A_optim_fr])^2
  
  VhatArray_A <- array(
        apply(V_sq_adj[1:num_y, 1:num_y, , drop = FALSE], 3,
              function(X) {
                X %*% as.matrix(A_means_used[1:num_y, 1:num_y]) %*% t(X)
                }
        ),
        c(num_y, num_y, NN_TT)
  )
  VhatArrayBdiagByTime <- IIGpkg::bdiagByTime(VhatArray_A, num_y, NN, TT, NN_num_y)

  out <- RcppSMCkalman::kfLGSSM(yObs = yObs,
                                uReg = NULL,
                                wReg = wReg,
                                A = Phi,
                                B = NULL,
                                C = B_means_used,
                                D = ccc,
                                Q = Q_used,
                                R = VhatArrayBdiagByTime,
                                initX = initX,
                                initU = initU,
                                initP = initP,
                                computeMFD = FALSE,
                                computeMSD = FALSE,
                                computePDD = TRUE,
                                computeLLH = TRUE)$kfLLHout
  if (out == Inf || is.na(out)) out <- -Inf
  if (PLOT) return(out)
  return(-out)
}
#' Helper to generate grid around a scalar parameter (centered at true value).
#'
#' @param true_val true parameter value
#' @param lower either `NULL` or lower bound
#' @param upper either `NULL` or upper bound
#' @param settings a named list of two elements, the first giving the number of
#'   grid points as `gridlength = 100` and the second the scaling `grid_width =
#'   0.9` which is used as follows: if both `lower` and `upper` are `NULL`, then
#'   the lower bound is computed as `true_val - true_val * settings$grid_width`
#'   and the upper bound is computed as
#'   `true_val - true_val * settings$grid_width`
#'
#' @return a numeric vector of grid points
#' @export
get_seq_grid_ml <- function(true_val, lower = NULL, upper = NULL,
                            settings = list(gridlength = 100,
                                            grid_width = 0.9)) {
  if (!is.null(lower)) {
    from_taken <- lower
  } else {
    from_taken <- true_val - true_val * settings$grid_width
  }
  if (!is.null(upper)) {
    to_taken <- upper
  } else {
    to_taken <- true_val + true_val * settings$grid_width
  }
  seq(from = from_taken,
      to = to_taken,
      length.out = settings$gridlength)
}
#' ML estimation under first identification scheme
#'
#' @param true_vals vector of true parameter values
#' @param start_vals starting values for optimization
#' @param HESSIAN logical; if `TRUE` then hessian is computed at maxima
#' @param arg_list arg list for all arguments (except `x`) of 
#'   [IIGpkg::ml_objective_identification_02]
#'
#' @return a list of two named elements, the results of the optimization as 
#'   `out_optim` and the print summary `summary_optim`
#' @export
ml_run_identification_01 <- function(true_vals, start_vals, HESSIAN = FALSE,
                                     arg_list = arg_list_used) {
  NN_num_y <- arg_list_used$NN * arg_list_used$num_y
  nidiofac <- NN_num_y
  num_a_optim <- length(arg_list$id_A_mat_optim)
  lower_bound <- c(0.0001,
                 rep(-Inf, times = NN_num_y - 1),
                 rep(0.0001, times = nidiofac),
                 rep(0.112, times = num_y),
                 rep(-Inf, times = num_a_optim - num_y))
  upper_bound <- rep(Inf, times = NN_num_y + nidiofac + num_a_optim)
  hessian_at_optim_vals <- NULL
  out_optim <- optimx::optimr(par = start_val,
                              fn = ml_objective_identification_01,
                              hess = NULL,
                              method = "L-BFGS-B",
                              lower = lower_bound,
                              upper = upper_bound,
                              hessian = HESSIAN,
                              control = list(trace = 3),
                              yObs = arg_list$yObs,
                              wReg = arg_list$wReg,
                              B_means_used = arg_list$B_used,
                              A_means_used = arg_list$A_used,
                              Q_used = arg_list$Q_used,
                              V_sq_adj = arg_list$V_sq_adj,
                              Phi = arg_list$Phi,
                              ccc = arg_list$ccc,
                              initX = arg_list$initX,
                              initU = arg_list$initU,
                              initP = arg_list$initP,
                              id_A_mat_optim = arg_list$id_A_mat_optim,
                              NN = arg_list$NN,
                              TT = arg_list$TT,
                              num_y = arg_list$num_y,
                              NN_TT = arg_list$NN_TT,
                              num_jnt_fac = arg_list$njfacs)
  optim_vals <- out_optim$par
  if (HESSIAN) hessian_at_optim_vals <- out_optim$hessian
  old <- options(pillar.sigfig = 6)
  # ci_scale <- sqrt(diag(solve(check_optim$hessian)))*1.96
  # ci <- cbind(check_optim$par - ci_scale,
  #             check_optim$par + ci_scale)
  # contained <- (true_vals >= ci[, 1] & true_vals <= ci[, 2])
  summary_optim <- tibble::tibble(`True value` = true_vals,
                                  `ML estiamte 01` = optim_vals,
                                  # `ML estiamte 02` = optim_vals2,
                                  `Start. val.` = start_vals) #,
                                  # `CI lower` = ci[, 1],
                                  # `CI upper` = ci[, 2],
                                  # contained)
  print(summary_optim, n = 999)
  options(old) 
  return(list(out_optim = list(out_optim = out_optim,
                               optim_vals = optim_vals,
                               hessian_at_optim_vals = hessian_at_optim_vals),
              summary_optim = summary_optim))
  #   

# check_optim$par
# check_optim2$par
# 
# # true_vals[(NN_num_y):(NN_num_y + num_y - 1)] <- A_means[1:num_y, 1:num_y][id_A_mat_optim[1:num_y]]
# optim_vals1 <- check_optim$par
# optim_vals2 <- check_optim2$par
# # ci_scale <- sqrt(diag(solve(check_optim$hessian)))*1.96
# # ci <- cbind(check_optim$par - ci_scale,
# #             check_optim$par + ci_scale)
# # contained <- (true_vals >= ci[, 1] & true_vals <= ci[, 2])
# old <- options(pillar.sigfig = 6)
# summary_optim <- tibble::tibble(`True value` = true_vals,
#                                 `ML estiamte 01` = optim_vals,
#                                 `ML estiamte 02` = optim_vals2,
#                                 `Start. val.` = start_val) #,
#                                 # `CI lower` = ci[, 1],
#                                 # `CI upper` = ci[, 2],
#                                 # contained)
# print(summary_optim, n = 999)
# options(old)
}
#' ML estimation under second identification scheme
#'
#' @param true_vals vector of true parameter values
#' @param start_vals starting values for optimization
#' @param HESSIAN logical; if `TRUE` then hessian is computed at maxima
#' @param arg_list arg list for all arguments (except `x`) of 
#'   [IIGpkg::ml_objective_identification_02]
#'
#' @return a list of two named elements, the results of the optimization as 
#'   `out_optim` and the print summary `summary_optim`
#' @export
ml_run_identification_02 <- function(true_vals, start_vals, HESSIAN = FALSE,
                                     arg_list = arg_list_used) {
  hessian_at_optim_vals <- NULL
  out_optim <- optim(par = start_val,
                     fn = ml_objective_identification_02,
                     gr = NULL,
                     method = "BFGS",
                     hessian = HESSIAN,
                     control = list(trace = 3),
                     yObs = arg_list$yObs,
                     wReg = arg_list$wReg,
                     B_means_used = arg_list$B_used,
                     A_means_used = arg_list$A_used,
                     Q_used = arg_list$Q_used,
                     V_sq_adj = arg_list$V_sq_adj,
                     Phi = arg_list$Phi,
                     ccc = arg_list$ccc,
                     initX = arg_list$initX,
                     initU = arg_list$initU,
                     initP = arg_list$initP,
                     id_A_mat_optim = arg_list$id_A_mat_optim,
                     NN = arg_list$NN,
                     TT = arg_list$TT,
                     num_y = arg_list$num_y,
                     NN_TT = arg_list$NN_TT,
                     num_jnt_fac = arg_list$njfacs
                     )
  optim_vals <- out_optim$par
  if (HESSIAN) hessian_at_optim_vals <- out_optim$hessian
  old <- options(pillar.sigfig = 6)
  # ci_scale <- sqrt(diag(solve(check_optim$hessian)))*1.96
  # ci <- cbind(check_optim$par - ci_scale,
  #             check_optim$par + ci_scale)
  # contained <- (true_vals >= ci[, 1] & true_vals <= ci[, 2])
  summary_optim <- tibble::tibble(`True value` = true_vals,
                                  `ML estiamte 01` = optim_vals,
                                  # `ML estiamte 02` = optim_vals2,
                                  `Start. val.` = start_vals) #,
                                  # `CI lower` = ci[, 1],
                                  # `CI upper` = ci[, 2],
                                  # contained)
  print(summary_optim, n = 999)
  options(old) 
  return(list(out_optim = list(out_optim = out_optim,
                               optim_vals = optim_vals,
                               hessian_at_optim_vals = hessian_at_optim_vals),
              summary_optim = summary_optim))
}