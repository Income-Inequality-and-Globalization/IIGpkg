#' Read and Preprocess Data
#'
#' Reads data from a specified path and filename, and preprocesses it by
#' removing specified variables.
#' 
#' @param out_gini_estimations an array of estimated Gini-coefficient values:
#'   rownmaes give the country (in the same order as the dataset), colnames the
#'   time, and the third dimension has the Bayesian a-posteriori mean, lower 
#'   and upper confidence band values
#' @param pth_base The base path to the data file.
#' @param fnm_data The name of the data file.
#' @param var_name_remove The name of the variable to be removed from the
#'   dataset.
#' @return A tibble with the specified variable removed.
#' @export
read_and_preprocess_data <- function(out_gini_estimations, 
                                     pth_base,
                                     fnm_data,
                                     var_name_remove) {
  data_raw <- read.csv(file.path(pth_base, fnm_data))
  data_raw <- data_raw %>% 
    tibble::as_tibble() %>%
    dplyr::select(-tidyselect::all_of(var_name_remove))

  gini_mean   <- as.vector(t(out_gini_estimations[, , "mean"])) * 100
  gini_ki_low <- as.vector(t(out_gini_estimations[, , "ki_low"])) * 100
  gini_ki_upp <- as.vector(t(out_gini_estimations[, , "ki_upp"])) * 100
  data_raw$gini_mean <- gini_mean
  data_raw$gini_ki_low <- gini_ki_low
  data_raw$gini_ki_upp <- gini_ki_upp
  return(data_raw)
}

#' Reshape Data to Long Format
#'
#' Reshapes the data into a long format suitable for ggplot2, focusing on specified Gini coefficients.
#' @param data_raw The raw data as a tibble.
#' @param ginis_to_use A character vector of the names of Gini coefficients to be used.
#' @return A tibble in long format with columns for Gini types and values.
#' @export
#' @examples
#' data_raw <- read_and_preprocess_data("/path/to/data", "data.csv", "variable_to_remove")
#' reshape_data_long(data_raw, c("gini1", "gini_from_gmm"))
reshape_data_long <- function(data_raw, ginis_to_use) {
  data_long <- data_raw %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(ginis_to_use), 
                 names_to = "gini_type", 
                 values_to = "gini_value")
  return(data_long)
}

#' Create Plot for a Given Country
#'
#' Creates a ggplot object for a given country, plotting Gini coefficient values over years.
#' @param data_long The data in long format.
#' @param name_country character giving the name of the country for which the 
#'    plot will be generated
#' @return A ggplot object.
#' @export
create_country_plot <- function(data_long, name_country) {
  data_long2 <- data_long[data_long$country == name_country, ]
  
  p <- ggplot2::ggplot(data_long2, ggplot2::aes(x = year)) +
    ggplot2::geom_point(ggplot2::aes(y = gini_value, color = gini_type)) + 
    ggplot2::geom_line(ggplot2::aes(y = gini_mean, color = "blue")) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = gini_ki_upp, ymin = gini_ki_low),
                         fill = "grey", alpha = 0.25) +
    ggplot2::labs(title = name_country,
                  y = "Gini value ( % )") +
    ggthemes::theme_tufte() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2")
  
  return(p)
}
