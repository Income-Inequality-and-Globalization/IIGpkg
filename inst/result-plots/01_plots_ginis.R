library(ggplot2)
library(dplyr)
library(tidyr)
library(IIGpkg)

# Main execution block
pth_dout <- file.path("./tests/testthat", "_data_test_measures")
test_dat <- readRDS(file.path(pth_dout, "output_gibbs_small.rds"))
# test_out <- subset_GibbsOutputIIG(test_dat$out_gibbs, 1:1000)
test_out <- readRDS("~/Dropbox/projects/IIG/IIGpkg/data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
test_out <- subset_GibbsOutputIIG(test_out, 100000:200000)
test_out <- compute_thin_mcmc(test_out, 10)
y_standardization  <- test_dat$yObs_standardization
y_centered_values  <- test_dat$yObs_centered_values
x_standardization  <- test_dat$x_standardization
x_centered_values  <- test_dat$x_centered_values
regs_to_use <- c("cpi_change", "unemployment", "gdp_ppp")
# regs_to_use <- c("gdp_ppp")
transformation_infos <- list(
  y_settings = list(
    centering = y_centered_values,
    scaling = y_standardization),
  x_settings = list(
    centering = x_centered_values,
    scaling = x_standardization,
    grid_length = 30,
    cut_off_num = 0,
    names_regs = c("cpi_change", "unemployment", "gdp_ppp")
  )
)
prc_measures <- precompute_measures_me(test_out, transformation_infos)
out_measures <- generate_measures_me(prc_measures$par_post_estim,
                                     prc_measures$post_me,
                                     regs_to_use,
                                     transformation_infos)
out_ginis_TT <- out_measures$gini_info_TT
out_ginis_KK <- out_measures$gini_info_KK
out_mus_TT   <- out_measures$mu_info_TT
out_mus_KK   <- out_measures$mu_info_KK
out_reg_grid <- prc_measures$regressor_grid_transformed
# create_me_plots_individual(out_ginis_KK, regs_to_use,
#                            settings = list(name_measure = "Gini",
#                                          	 plot_grid = c(4, 5),
#                                          	 WITH_CI = TRUE))
# testme <- create_me_plots_individual(out_mus_KK, regs_to_use,
#                            settings = list(name_measure = "Expt_Inc_E[Y]",
#                                            plot_grid = c(4, 5),
#                                            plot_type = "ggplot",
#                                            WITH_CI = TRUE))
out_reg_grid <- prc_measures$regressor_grid_transformed

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
                      ADD_KI = TRUE,
                      grid_dim = c(4, 5))
out_measure_gini <- list(info_TT = out_ginis_TT, info_KK = out_ginis_KK)
plot_ginis_TT <- generate_country_plots(pth_data,
                                        output_measure = out_measure_gini,
                                        output_regs_grid = out_reg_grid,
                                        vars_gini,
                                        settings_gini,
                                        "Gini")
vars_mu <- "est_mean"
settings_mu <- list(scale_measure = 100,
                    y_lab = "Expt. Inc. E[Y]",
                    x_lab = "time",
                    x_var = NULL,
                    ADD_KI = TRUE,
                    grid_dim = c(4, 5))
out_measure_mu <- list(info_TT = out_mus_TT, info_KK = out_mus_KK)
plot_mus_TT <- generate_country_plots(pth_data,
                                      output_measure = out_measure_mu,
                                      output_regs_grid = out_reg_grid,
                                      vars_mu,
                                      settings_mu,
