library(ggplot2)
library(dplyr)
library(tidyr)
library(IIGpkg)

# Main execution block
pth_dout <- file.path("./tests/testthat", "_data_test_measures")
test_dat <- readRDS(file.path(pth_dout, "output_gibbs_small.rds"))
# test_out <- subset_GibbsOutputIIG(test_dat$out_gibbs, 1:1000)
# test_out <- readRDS("/home/iz/Documents/save_good_model/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj0_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj0_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")

# Models without hyperpriors
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi1000_nu1000_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi1000_nu1000_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7_IO5000.rds")

# 1000-1000 type models in three versions: 5/6, 10/3 and 10/0.1
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi1000_nu1000/mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi1000_nu1000mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi1000_nu1000/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi1000_nu1000mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi1000_nu1000/mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi1000_nu1000mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi1000_nu1000/mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi1000_nu1000mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000.rds")

# 27-7 type models in three versions: 5/6, 10/3 and 10/0.1
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi27_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000_backup/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi27_nu7mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000.rds")

# 57-7 type models in three versions: 5/6, 10/3 and 10/0.1
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7mu_b0_0_Sigma_b0_100_alpha_b0_255_beta_b0_256_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7mu_b0_0_Sigma_b0_100_alpha_b0_5_beta_b0_6_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7mu_b0_0_Sigma_b0_100_alpha_b0_3_beta_b0_10_IO5000.rds")
# test_out <- readRDS("./data/output/empirical/Estimates/cs1_pj1_Reg7_B0_Om1000_D0_OmD1000_A1_Psi57_nu7/mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000/cs1_pj1_B0_Om_D0_OmD1000_A1_Psi57_nu7mu_b0_0_Sigma_b0_100_alpha_b0_10_beta_b0_0.1_IO5000.rds")

test_out <- subset_GibbsOutputIIG(test_out, 100000:200000)
test_out <- compute_thin_mcmc(test_out, 10)
y_standardization  <- test_dat$yObs_standardization
y_centered_values  <- test_dat$yObs_centered_values
x_standardization  <- test_dat$x_standardization
x_centered_values  <- test_dat$x_centered_values
wRegRawList        <- test_dat$wRegRawList

testme <- read.table("./data/input/data-sm/data_coef_covariates.txt", header = TRUE)
testme$gdp_ppp_pc[testme$country != "Venezuela"] <- testme$gdp_ppp[testme$country != "Venezuela"] / testme$population[testme$country != "Venezuela"]
testme$gdp_ppp_pc[testme$country == "Venezuela"] <- testme$gdp_ppp[testme$country == "Venezuela"]
testme$gdp_ppp_pc[testme$country == "Venezuela"] <- testme$gdp_ppp[testme$country == "Venezuela"]
testme2 <- testme %>% dplyr::filter(!year %in% 2021)
unq_cont <- unique(testme$country)
id_repl <- seq(from = 3, by = 3, length.out = 10)
for (i in 1:length(unq_cont)) {
  wRegRawList[id_repl[i], ] <- testme2  %>% dplyr::filter(country == unq_cont[i]) %>% dplyr::pull(gdp_ppp_pc) / 1000
}
regs_to_use <- c("cpi_change", "unemployment", "gdp_ppp")
# regs_to_use <- c("gdp_ppp")
transformation_infos <- list(
  y_settings = list(
    centering = y_centered_values,
    scaling = y_standardization),
  x_settings = list(
    grid_length = 30,
    cut_off_num = 0,
    names_regs = c("cpi_change", "unemployment", "gdp_ppp")
  ),
  meta = list(
    wRegs_raw = wRegRawList,
    wRegs_grd = test_out$initials$wReg
  )
)
prc_measures <- precompute_measures_me(test_out, transformation_infos)
out_measures <- generate_measures_me(
  prc_measures$post_par_estimates,
  prc_measures$post_par_estimates_forward_diff,
  prc_measures$post_me_grid,
  regs_to_use,
  transformation_infos)

out_ginis_TT    <- out_measures$gini_info_TT
out_ginis_KK    <- out_measures$gini_info_KK
out_ginis_ft_TT <- out_measures$gini_info_fd_TT

out_mus_TT    <- out_measures$mu_info_TT
out_mus_KK    <- out_measures$mu_info_KK
out_mus_ft_TT <- out_measures$mu_info_fd_TT

out_reg_grid <- list(
  grid_KK = prc_measures$regressor_grid$reg_grid_original,
  # grid_KK = prc_measures$regressor_grid$reg_grid_transformed)
  grid_TT = transformation_infos$meta$wRegs_raw)

# create_me_plots_individual(out_ginis_KK, regs_to_use,
#                            settings = list(name_measure = "Gini",
#                                          	 plot_grid = c(4, 5),
#                                            plot_type = "base",
#                                          	 WITH_CI = TRUE))
# testme <- create_me_plots_individual(out_mus_KK, regs_to_use,
#                            settings = list(name_measure = "Expt_Inc_E[Y]",
#                                            plot_grid = c(4, 5),
#                                            plot_type = "base",
#                                            WITH_CI = TRUE))

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
settings_gini_me <- list(scale_measure = 1,
                         y_lab = "ME on Gini",
                         x_lab = "time",
                         x_var = NULL,
                         ADD_KI = TRUE)
settings_gini = list(
  measure = settings_gini,
  measure_me = settings_gini_me
)
out_measure_gini <- list(info_TT = out_ginis_TT,
                         info_KK = out_ginis_KK,
                         info_FD = out_ginis_ft_TT)
plot_ginis_TT <- generate_country_plots(pth_data,
                                        output_measure = out_measure_gini,
                                        output_regs_grid = out_reg_grid,
                                        vars_gini,
                                        regs_to_use,
                                        settings_gini,
                                        "Gini",
                                        grid_dim = c(4, 5))
vars_mu <- "est_mean"
settings_mu <- list(scale_measure = 1,
                    y_lab = "Expt. Inc. E[Y]",
                    x_lab = "time",
                    x_var = NULL,
                    ADD_KI = TRUE)
settings_mu_me <- list(scale_measure = 1,
                       y_lab = "ME on E[Y]",
                       x_lab = "time",
                       x_var = NULL,
                       ADD_KI = TRUE)
out_measure_mu <- list(info_TT = out_mus_TT,
                       info_KK = out_mus_KK,
                       info_FD = out_mus_ft_TT)
settings_mu = list(
  measure = settings_mu,
  measure_me = settings_mu_me
)
plot_mus_TT <- generate_country_plots(pth_data,
                                      output_measure = out_measure_mu,
                                      output_regs_grid = out_reg_grid,
                                      vars_mu,
                                      regs_to_use,
                                      settings_mu,
                                      "Exp.Inc",
                                      grid_dim = c(4,5))
