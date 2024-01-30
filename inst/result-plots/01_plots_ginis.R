library(ggplot2)
library(dplyr)
library(tidyr)
library(IIGpkg)

# Main execution block
pth_dta <- file.path("./tests/testthat", "_data_test_measures")
test_data <- readRDS(file.path(pth_dta, "output_gibbs_small.rds"))
test_out <- test_data$out_gibbs
y_standardization  <- test_data$yObs_standardization
y_centered_values  <- test_data$yObs_centered_values
transformation_infos <- list(centering = y_centered_values,
                             scaling = y_standardization,
                             names_regs = c("cpi_change",
                                            "unemployment",
                                            "gdp_ppp"))
out_measures <- generate_measures_me(test_out, transformation_infos)
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