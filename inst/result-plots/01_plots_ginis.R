library(ggplot2)
library(dplyr)
library(tidyr)
library(IIGpkg)

# Main execution block
pth_dta <- file.path("./tests/testthat", "_data_test_measures")
test_data <- readRDS(file.path(pth_dta, "output_gibbs_small.rds"))
test_out  <- test_data$out_gibbs
y_standardization  <- test_data$yObs_standardization
y_centered_values  <- test_data$yObs_centered_values
transformation_infos <- list(
  y_settings = list(
    centering = y_centered_values,
    scaling = y_standardization),
  x_settings = list(
    grid_length = 100,
    cut_off_num = 4,
    names_regs = c("cpi_change",
                   "unemployment",
                   "gdp_ppp")
  )
)
out_measures <- generate_measures_me(test_out, transformation_infos)

# out_measures$gini_info_KK$cpi_change <- out_measures$gini_info_KK$cpi_change[, , 2, ]
# out_measures$gini_info_KK$unemployment <- out_measures$gini_info_KK$unemployment[, , 2, ]
# out_measures$gini_info_KK$gdp_ppp <- out_measures$gini_info_KK$gdp_ppp[, , 2, ]
# out_measures$regressor_grid_transformed <- out_measures$regressor_grid_transformed[,, 2,]

var_tkn <- "unemployment"
par(mfrow = c(3, 3))
nn_tk <- 1
for (nn_tk in 1:9) {
    # tmp_test <- apply(out_measures$gini_info_KK$unemployment[, , ,], c(1, 3, 4), mean)
    #   testme_mean <- tmp_test[nn_tk, , 1]
    # testme_kiup <- tmp_test[nn_tk, , 2]
    # testme_kilo <- tmp_test[nn_tk, , 3]
    # 
    # min_y <- min(tmp_test[nn_tk, , ])
    # max_y <- max(tmp_test[nn_tk, , ])
    # plot(testme_mean, type = "l", ylim = c(min_y, max_y))
    # lines(testme_kiup, col = "blue")
    # lines(testme_kilo, col = "blue")
  tmp_test <- out_measures$gini_info_KK[[var_tkn]][, , ,]
  min_y <- min(tmp_test[nn_tk, , , 1])
  max_y <- max(tmp_test[nn_tk, , , 1])
  plot(tmp_test[nn_tk, 1, , 1], type = "l", ylim = c(min_y, max_y))
  for (tt in 2:21) {
    testme_mean <- tmp_test[nn_tk, tt, ,1]
    lines(testme_mean)
    # testme_kiup <- tmp_test[nn_tk, tt, ,2]
    # testme_kilo <- tmp_test[nn_tk, tt, ,3]
    # 
    # 
    # plot(testme_mean, type = "l", ylim = c(min_y, max_y))
    # lines(testme_kiup, col = "b")
    # lines(testme_kilo, col = "blue")
  }
}

par(mfrow = c(3, 3))
nn_tk <- 1
for (nn_tk in 1:9) {
  tmp_test <- apply(out_measures$gini_info_KK[[var_tkn]][, , ,], c(1, 3, 4), mean)
    testme_mean <- tmp_test[nn_tk, , 1]
  testme_kiup <- tmp_test[nn_tk, , 2]
  testme_kilo <- tmp_test[nn_tk, , 3]

  min_y <- min(tmp_test[nn_tk, , ])
  max_y <- max(tmp_test[nn_tk, , ])
  plot(testme_mean, type = "l") #, ylim = c(min_y, max_y))
  # lines(testme_kiup, col = "blue"#)
  #lines(testme_kilo, col = "blue")
}

# min_y <- min(out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , ])
# max_y <- max(out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , ])
plot(testme_mean, type = "l") #, ylim = c(min_y, max_y))

nn_tk <- 1
par(mfrow = c(4, 5))
for (tt_tk in 1:20) {
  testme_mean <- out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , 1]
  testme_kiup <- out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , 2]
  testme_kilo <- out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , 3]
  
  # min_y <- min(out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , ])
  # max_y <- max(out_measures$gini_info_KK$gdp_ppp[nn_tk, tt_tk, , ])
  plot(testme_mean, type = "l") #, ylim = c(min_y, max_y))
  # lines(testme_kiup, col = "blue")
  # lines(testme_kilo, col = "blue")
}


out_ginis_TT <- out_measures$gini_info_TT
out_ginis_KK <- out_measures$gini_info_KK
out_mus_TT   <- out_measures$mu_info_TT
out_mus_KK   <- out_measures$mu_info_KK
# out_reg_grid <- out_measures$regressor_gird_original
out_reg_grid <- out_measures$regressor_grid_transformed

pth_base <- "/home/iz/Dropbox/projects/IIG/IIGpkg/data/input/data-sm"
pth_data <- file.path(pth_base, "data_coef_covariates_with_ginis_2020.csv")

names_measures_gini <- list(mean = "gini_mean",
                            ki_lower = "gini_ki_low",
                            ki_upper = "gini_ki_upp")
vars_gini <- setdiff(c("gini1", "gini2", "gini_from_gmm"), c("gini1", "gini2"))
settings_gini <- list(scale_measure = 100,
                      y_lab = "Gini ( % )",
                      x_lab = "time",
                      x_var = NULL,
                      ADD_KI = TRUE)
plot_ginis_TT <- generate_country_plots(pth_data,
                                        out_ginis_TT,
                                        NULL,
                                        vars_gini,
                                        settings_gini,
                                        "gini")
settings_gini$ADD_KI <- FALSE
settings_gini$X_TRANSFORMED <- TRUE

settings_gini$x_lab <- "CPI change"
settings_gini$x_var <- "cpi_change"
plot_ginis_cpi_change <- generate_country_plots(pth_data,
                                                out_ginis_KK[["cpi_change"]],
                                                out_reg_grid,
                                                vars_gini,
                                                settings_gini,
                                                "gini")
settings_gini$x_lab <- "Unemployment"
settings_gini$x_var <- "unemployment"
plot_ginis_unemployment <- generate_country_plots(pth_data,
                                                  out_ginis_KK[["unemployment"]],
                                                  out_reg_grid,
                                                  vars_gini,
                                                  settings_gini,
                                                  "gini")
settings_gini$x_lab <- "GDP"
settings_gini$x_var <- "gdp_ppp"
plot_ginis_gdp_ppp <- generate_country_plots(pth_data,
                                             out_ginis_KK[["gdp_ppp"]],
                                             out_reg_grid,
                                             vars_gini,
                                             settings_gini,
                                             "gini")


mus_to_use <- "est_mean"
settings_mu <- list(scale_measure = 1, y_lab = "Expt. Income ( E[Y] )")
plot_mus <- generate_country_plots(pth_data,
                                   out_mus_TT,
                                   mus_to_use,
                                   settings_mu,
                                   "mu")