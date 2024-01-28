library(ggplot2)
library(dplyr)
library(tidyr)
library(IIGpkg)

# Reshaping the data to long format for ggplot

# Main execution block
pth_dta <- file.path("./tests/testthat",
                     "_data_test_measures")
test_data <- readRDS(file.path(pth_dta, "output_gibbs_small.rds"))
test_out <- test_data$out_gibbs
y_standardization  <- test_data$yObs_standardization
y_centered_values  <- test_data$yObs_centered_values
transformation_infos <- list(centering = y_centered_values,
                             scaling = y_standardization)
out_measures <- generate_measures_me(test_out, transformation_infos)
out_ginis <- out_measures$gini_info
out_mus <- out_measures$mu_info



pth_base <- "/home/iz/Dropbox/projects/IIG/IIGpkg/data/input/data-sm"
pth_data <- file.path(pth_base, 
                      "data_coef_covariates_with_ginis_2020.csv")

var_name_remove <- c("gini1", "gini2")
data_ginis <- data_plot_processing(out_ginis, 
                                   pth_data,
                                  list(mean = "gini_mean",
                                       ki_lower = "gini_ki_low",
                                       ki_upper = "gini_ki_upp"),
                                  var_name_remove,
                                  100)
ginis_to_use <- setdiff(c("gini1", "gini2", "gini_from_gmm"), var_name_remove)
data_long_ginis <- get_data_plot_long(data_ginis,
                                      ginis_to_use,
                                      "gini_type",
                                      "gini_value")


data_mus <- data_plot_processing(out_mus, 
                                 pth_data,
                                 list(mean = "mu_mean",
                                      ki_lower = "mu_ki_low",
                                      ki_upper = "mu_ki_upp"),
                                 NULL)
mus_to_use <- "est_mean"
data_long_mus <- get_data_plot_long(data_mus,
                                    mus_to_use,
                                    "mu_type",
                                    "mu_value")
unique_countries <- unique(data_long_mus$country)
list_plots_per_country <- list()
for (country in unique_countries[-c(length(unique_countries))]) {
  list_plots_per_country[[country]] <- create_single_country_plot(
    data_long = data_long_mus,
    data_info = list(
      name_country = country,
      y_lab = "E[Y] = mu of SM-Income-Distr."
    ),
    measure_info = list(
      vals = "mu_value",
      types = "mu_type"
    ),
    estim_infos = list(
      col_mean = "mu_mean",
      col_ki_upp = "mu_ki_upp",
      col_ki_low = "mu_ki_low"
    )
  )
}
gridExtra::grid.arrange(grobs = list_plots_per_country, as.table = FALSE)
my_grobs <- lapply(list_plots_per_country, ggplot2::ggplotGrob)

unique_countries <- unique(data_long_ginis$country)
list_plots_per_country <- list()
for (country in unique_countries[-c(length(unique_countries))]) {
  list_plots_per_country[[country]] <- create_single_country_plot(
    data_long = data_long_ginis,
    data_info = list(
      name_country = country,
      y_lab = "Gini ( % )"
    ),
    measure_info = list(
      vals = "gini_value",
      types = "gini_type"
    ),
    estim_infos = list(
      col_mean = "gini_mean",
      col_ki_upp = "gini_ki_upp",
      col_ki_low = "gini_ki_low"
    )
  )
}
gridExtra::grid.arrange(grobs = list_plots_per_country, as.table = FALSE)
my_grobs <- lapply(list_plots_per_country, ggplot2::ggplotGrob)