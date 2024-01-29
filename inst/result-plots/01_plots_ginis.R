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

pth_base <- "/home/iz/Dropbox/projects/IIG/IIGpkg/data/input/data-sm"
pth_data <- file.path(pth_base, "data_coef_covariates_with_ginis_2020.csv")

names_measures_gini <- list(mean = "gini_mean",
                            ki_lower = "gini_ki_low",
                            ki_upper = "gini_ki_upp")
vars_gini <- setdiff(c("gini1", "gini2", "gini_from_gmm"), c("gini1", "gini2"))
settings_gini <- list(scale_transform = 100,
                      y_lab = "Gini ( % )",
                      x_lab = "time")
plot_ginis_TT <- generate_country_plots(pth_data,
                                        out_ginis_TT,
                                        vars_gini,
                                        settings_gini,
                                        "gini")
settings_gini$x_lab <- "CPI change"
plot_ginis_cpi_change <- generate_country_plots(pth_data,
                                                out_ginis_KK[["cpi_change"]],
                                                vars_gini,
                                                settings_gini,
                                                "gini")
settings_gini$x_lab <- "Unemployment"
plot_ginis_unemployment <- generate_country_plots(pth_data,
                                                  out_ginis_KK[["unemployment"]],
                                                  vars_gini,
                                                  settings_gini,
                                                  "gini")
settings_gini$x_lab <- "GDP"
plot_ginis_gdp_ppp <- generate_country_plots(pth_data,
                                             out_ginis_KK[["gdp_ppp"]],
                                             vars_gini,
                                             settings_gini,
                                             "gini")


mus_to_use <- "est_mean"
settings_mu <- list(scale_transform = 1, y_lab = "Expt. Income ( E[Y] )")
plot_mus <- generate_country_plots(pth_data,
                                    out_mus_TT,
                                    mus_to_use,
                                    settings_mu,
                                    "mu")