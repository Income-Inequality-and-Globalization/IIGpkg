#' Plotting function for fitted values and their components.
#' 
#' Plots observations (data/measurements), the fit coming from joint factors and
#' loading, the fit stemming from idiosyncratic factors and loadings, the fit
#' stemming from regressors (times their marginal effects), and the overall fit
#' to the data (evaluated at the posterior means).
#'
#' @param yObs observations (data/measurements)
#' @param wRegs regressors
#' @param f_pm fitted factor (a posteriori) means; possibly after burnin
#' @param B_pm fitted loading (a posteriori) means; possibly after burnin
#' @param D_pm fitted regressor marginal (a posteriori) means; possibly after burnin
#' @param NN_num_y number of cross sectional units times dimension of the data
#' @param settings a named list of settings; first element is the path where 
#'   plots are saved and the second is the grid according to which plots are 
#'   arranged i.e. `list(path_save = ".", plot_grid = c(3, 3)`
#'
#' @return pure side effect function returning invisibly `NULL`
#' @export
plot_fit_diagnostics <- function(yObs, wRegs, f_pm, B_pm, D_pm, NN_num_y, 
                                 settings = list(path_save = ".",
                                                 plot_grid = c(3, 3))) {
  plot_title <- paste0("black: data;   ",
                       "purple: joint fit;   ",
                       "red: B*f_jnt:   ",
                       "orange: B*f_idio;   ",
                       "green: D*w_regs;   ")
  iter_plot <- 1
  break_plot_num <- prod(settings$plot_grid)
  pn <- get_plot_name(settings$path_save, "fitted_plot_num", iter_plot)
  par(mfrow = settings$plot_grid)
  for (ii in seq_len(NN_num_y)) {
    if (ii %% break_plot_num == 1) {pdf(pn);par(mfrow = settings$plot_grid)}
    f_Bjoint_fit <- f_pm[1, ] * B_pm[ii, 1]
    f_Bidios_fit <- f_pm[ii + 1, ] * B_pm[ii, ii + 1]
    D_wRegs_fit  <- (D_pm %*% w_regs)[ii, ]
    y_means_fit  <- f_Bjoint_fit + f_Bidios_fit + D_wRegs_fit

    lims <- c(yObs[ii, ], f_Bjoint_fit, f_Bidios_fit, D_wRegs_fit, y_means_fit)
    ylim_min <- min(lims, na.rm = TRUE)
    ylim_max <- max(lims, na.rm = TRUE)
  
    plot(yObs[ii, ], type = "l", ylim = c(ylim_min, ylim_max), lwd = 2,
         ylab =  paste0("No. ", ii, " time series"), xlab = "time")

    lines(f_Bjoint_fit, col = "red")
    lines(f_Bidios_fit, col = "orange")
    lines(D_wRegs_fit, col = "green")

    lines(y_means_fit, col = "purple")
    if (ii %% break_plot_num == 0) {
      mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
      iter_plot <- iter_plot + 1
      dev.off()
      pn <- get_plot_name(settings$path_save, "fitted_plot_num", iter_plot)
    }
  }
  mtext(plot_title, side = 3, line = -2, outer = TRUE)
  dev.off()
  return(invisible(NULL))
}
#' Plotting function for fitted vs. true factors.
#' 
#' Plots posteriori meaned factors vs. their true values.
#'
#' @param f_true true posterior factors
#' @inheritParams plot_fit_diagnostics
#' @inheritParams plot_fit_diagnostics
#' @param num_joint_fac 
#'
#' @return pure side effect function returning invisibly `NULL`
#' @export
plot_factor_diagnostics <- function(f_true, f_pm, NN_num_y, num_joint_fac,
                                    settings = list(path_save = ".",
                                                    plot_grid = c(3, 3))) {
  plot_title <- paste0("black: true simulated factor;   ",
                       "red: posterior factor;   ")
  iter_plot <- 1
  break_plot_num <- prod(settings$plot_grid)
  pn <- get_plot_name(settings$path_save, "factor_plot_num", iter_plot)
  par(mfrow = settings$plot_grid)
  for (ii in seq_len(NN_num_y + num_joint_fac)) {
    if (ii %% break_plot_num == 1) {pdf(pn);par(mfrow = settings$plot_grid)}

    ylim_min <- min(c(f_true[ii, ], f_pm[ii, ]), na.rm = TRUE)
    ylim_max <- max(c(f_true[ii, ], f_pm[ii, ]), na.rm = TRUE)
  
    plot(f_true[ii, ], type = "l", ylim = c(ylim_min, ylim_max),
         lwd = 2, ylab =  paste0("No. ", ii, " time series"), xlab = "time")
    lines(f_pm[ii, ], col = "red")
    if (ii %% break_plot_num == 0) {
      mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
      iter_plot <- iter_plot + 1
      dev.off()
      pn <- get_plot_name(settings$path_save, "factor_plot_num", iter_plot)
    }
  }
  mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
  dev.off()
  return(invisible(NULL))
}
#' Plotting function for B*f for idiosyncratic factors
#' 
#' Idiosyncratic factors: fitted vs. true multiplied by the loadings.
#'
#' @param f_true true posterior factors
#' @param B_true true loadings
#' @inheritParams plot_fit_diagnostics
#' @inheritParams plot_fit_diagnostics
#'
#' @return pure side effect function returning invisibly `NULL`
#' @export
plot_fBido_diagnostics <- function(f_true, f_pm, B_true, B_pm, NN_num_y,
                                   settings = list(path_save = ".",
                                                   plot_grid = c(3, 3))) {
  plot_title <- paste0("black: B*f true simulated idio;   ",
                       "green: B*f posterior factor idio;   ")
  iter_plot <- 1
  break_plot_num <- prod(settings$plot_grid)
  pn <- get_plot_name(settings$path_save, "Bf_idio_plot_num", iter_plot)
  par(mfrow = settings$plot_grid)
  for (ii in seq_len(NN_num_y)) {
    if (ii %% break_plot_num == 1) {pdf(pn);par(mfrow = settings$plot_grid)}
    fBtrue <- f_true[ii + 1, ] * B_true[ii, ii + 1]
    fBmean <- f_pm[ii + 1, ] * B_pm[ii, ii + 1]
    lims <- c(fBtrue, fBmean)
    ylim_min <- min(lims)
    ylim_max <- max(lims)
    plot(fBtrue, col = "black",
         type = "l", ylim = c(ylim_min, ylim_max),
         ylab = paste0("No. ", ii, " time series"))
    lines(fBmean, col = "green")

    if (ii %% break_plot_num == 0) {
      mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
      iter_plot <- iter_plot + 1
      dev.off()
      pn <- get_plot_name(settings$path_save, "Bf_idio_plot_num", iter_plot)
    }
  }
  mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
  dev.off()
  return(invisible(NULL))
}
#' Plotting function for B*f for idiosyncratic factors
#' 
#' Idiosyncratic factors: fitted vs. true multiplied by the loadings.
#'
#' @param f_true true posterior factors
#' @param B_true true loadings
#' @inheritParams plot_fit_diagnostics
#' @inheritParams plot_fit_diagnostics
#'
#' @return pure side effect function returning invisibly `NULL`
#' @export
plot_fBjnt_diagnostics <- function(f_true, f_pm, B_true, B_pm,
                                   NN_num_y, num_joint_fac,
                                   settings = list(path_save = ".",
                                                   plot_grid = c(3, 3))) {
  plot_title <- paste0("black: B*f true simulated idio;   ",
                       "green: B*f posterior factor idio;   ")
  iter_plot <- 1
  break_plot_num <- prod(settings$plot_grid)
  pn <- get_plot_name(settings$path_save, "Bf_common_plot_num", iter_plot)
  par(mfrow = settings$plot_grid)
  for (ii in seq_len(NN_num_y)) {
    if (ii %% break_plot_num == 1) {pdf(pn);par(mfrow = settings$plot_grid)}
    fBtrue <- f_true[1, ] * B_true[ii, 1:num_joint_fac]
    fBmean <- f_pm[1, ] * B_pm[ii, 1:num_joint_fac]
    lims <- c(fBtrue, fBmean)
    ylim_min <- min(lims)
    ylim_max <- max(lims)
    plot(fBtrue, col = "black",
         type = "l", ylim = c(ylim_min, ylim_max),
         ylab = paste0("No. ", ii, " time series"))
    lines(fBmean, col = "red")

    if (ii %% break_plot_num == 0) {
      mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
      iter_plot <- iter_plot + 1
      dev.off()
      pn <- get_plot_name(settings$path_save, "Bf_common_plot_num", iter_plot)
    }
  }
  mtext(plot_title, side = 3, line = -1.5, outer = TRUE)
  dev.off()
  return(invisible(NULL))
}
get_plot_name <- function(base_path, type_name, iter) {
  paste0(base_path, "/", type_name, "_",
         formatC(iter, digits = 1, format = "d", flag = "0"),
         ".pdf")
}


# plot_title <- paste0("black: B*f true simulated common;   ",
#                      "red: B*f posterior factor common;   ")
# B_times_f_true <- B_means %*% f_means_simul
# B_times_f_post <- B_means_sim_out %*% f_all
# iter_plot <- 1
# plot_name <- paste0("data/output/tmp", "Bf_common_plot_num_", iter_plot, ".pdf")
# for (t_taken in 1:30) {
#   if (t_taken %% 9 == 1) {
#     pdf(plot_name)
#     par(mfrow = c(3, 3))
#   }
#     ylim_min <- min(c(f_means_simul[1, ] * B_means[t_taken, 1], f_all[1, ] * B_means_sim_out[t_taken, 1]))
#     ylim_max <- max(c(f_means_simul[1, ] * B_means[t_taken, 1], f_all[1, ] * B_means_sim_out[t_taken, 1]))
#     plot(f_means_simul[1, ] * B_means[t_taken, 1], col = "black", type = "l",
#          ylim = c(ylim_min, ylim_max), ylab =  paste0("No. ", t_taken, " time series"))
#     lines(f_all[1, ] * B_means_sim_out[t_taken, 1], col = "red",
#           ylim = c(ylim_min, ylim_max))
#   if (t_taken %% 9 == 0) {
#     mtext(plot_title, side = 3, line = -2, outer = TRUE)
#     iter_plot <- iter_plot + 1
#     dev.off()
#     plot_name <- paste0("data/output/tmp", "Bf_common_plot_num_", iter_plot, ".pdf")
#   }
# }
# mtext(plot_title, side = 3, line = -2, outer = TRUE)
# dev.off()