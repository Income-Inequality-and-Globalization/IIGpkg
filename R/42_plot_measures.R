#' Read and Preprocess Data
#'
#' Reads data from a specified `path` and `filename`, and preprocesses it by
#' removing specified variables, and adding output measures from a suitable 
#' array (`out_measures`).
#' 
#' @param out_measures an array of estimated measure-coefficient values: rownmaes
#'   give the country (in the same order as the dataset), colnames the time, and
#'   the third dimension has the Bayesian a-posteriori mean, lower and upper
#'   confidence band values. The measure values can be Ginis, or the expectaion
#'   as estimated from the Bayesian parameter samples a-posteriori.
#' @param pth_data The path to the data file which appended by content of the
#'   of the first argument (ie. transform the arrary to the data format and 
#'   add its content as new variables)
#' @param names_measures A named list or vector specifying the column names to be 
#'   created in the dataset for each type of measure. It must include names for 
#'   "mean", "ki_low", and "ki_upp", corresponding to the Bayesian a-posteriori 
#'   mean, lower confidence band, and upper confidence band values, respectively. 
#'   The names are used to dynamically create and assign these measures to the 
#'   appropriate columns in the dataset. For example, names_measures could be 
#'   set as `list(mean = "gini_mean", ki_low = "gini_ki_low", 
#'   ki_upp = "gini_ki_upp")`.
#' @param var_name_remove The name of the variables to be removed from the
#'   dataset, if necessary.
#' @param scale_values A numeric value used to scale the measure values extracted
#'   from `out_measures`. This scaling factor is applied to the mean, lower, and 
#'   upper confidence band values. It allows for adjusting the measure values 
#'   according to the desired scale (e.g., converting proportions to percentages).
#'   Default is 1 (no scaling).
#'
#' @return A tibble with the specified variable removed, and the output measure
#'   variables added
#'
#' @export
data_plot_processing <- function(out_measures,
                                 pth_data,
                                 names_measures,
                                 names_var_remove = NULL,
                                 scale_values = 1) {
  data_raw <- read.csv(file.path(pth_data)) %>% tibble::as_tibble() 

  if (!is.null(var_name_remove)) {
    data_raw <- data_raw %>% 
      dplyr::select(-tidyselect::all_of(names_var_remove))
  }

  measure_mean_01 <- get_measure_vals(out_measures, "mean", scale_values)
  measure_ki_l_02 <- get_measure_vals(out_measures, "ki_low", scale_values)
  measure_ki_u_03 <- get_measure_vals(out_measures, "ki_upp", scale_values)
  data_raw[names_measures$mean]   <- measure_mean_01
  data_raw[names_measures$ki_low] <- measure_ki_l_02
  data_raw[names_measures$ki_upp] <- measure_ki_u_03
  return(data_raw)
}
get_measure_vals <- function(out_msr, name_var, scl_val = 1) {
  as.vector(t(out_msr[, , name_var])) * scl_val
}

#' Reshape Data to Long Format
#'
#' Reshapes the data into a long format suitable for ggplot2, focusing on
#' specified variable(s).
#'
#' This function transforms the data into a long format, which is particularly
#' useful for ggplot2 visualizations. It allows for the selection of specific
#' variables to be reshaped, making it versatile for various types of data beyond
#' just Gini coefficients.
#'
#' @param data_raw The raw data as a tibble.
#' @param var_to_use A character vector of variable names to be used in the
#'   transformation. These variables will be reshaped from wide to long format.
#' @param nms_to A string specifying the name of the new column in the long
#'   format that will contain the variable names from `var_to_use`.
#' @param vals_to A string specifying the name of the new column in the long
#'   format that will contain the values from the variables in `var_to_use`.
#' @return A tibble in long format with columns as specified by `nms_to` and
#'   `vals_to`.
#'
#' @export
#' @examples
#' \dontrun{
#' # Example usage with Gini coefficients
#' data_long_ginis <- get_data_gini_long(data_raw,
#'                                       ginis_to_use,
#'                                       "gini_type",
#'                                       "gini_value")
#' }
get_data_plot_long <- function(data_raw, var_to_use, nms_to, vals_to) {
  data_long <- data_raw %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(var_to_use), 
                        names_to = nms_to, 
                        values_to = vals_to)
  return(data_long)
}


#' Reshape Data to Long Format
#'
#' Reshapes the data into a long format suitable for ggplot2, focusing on
#' specified expectation variable.
#'
#' @param data_raw The raw data as a tibble.
#' @param ginis_to_use character giving the name of the `mu`/expectation
#'   variable to use
#' @return A tibble in long format with columns for Gini types and values.
#' @export
#' @examples
#' data_raw <- read_and_preprocess_data("/path/to/data", "data.csv", "variable_to_remove")
#' reshape_data_long(data_raw, c("gini1", "gini_from_gmm"))
get_data_mu_long <- function(data_raw, mu_to_use) {
  browser()
  data_long <- data_raw %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(mu_to_use), 
      names_to = "gini_type", 
      values_to = "mu_value")
  return(data_long)
}
#' Create Plot for a Given Country
#'
#' Creates a ggplot object for a given country, plotting specified measure
#' values over years.
#' 
#' @param data_long The data in long format.
#' @param data_info A list containing information about the data and plot:
#'   \itemize{
#'     \item \code{name_country}: character giving the name of the country for
#'     which the plot will be generated.
#'     \item \code{y_lab}: label for the y-axis.
#'   }
#' @param measure_info A list containing information about the measure points:
#'   \itemize{
#'     \item \code{vals}: The name(s) of the column(s) in \code{data_long}
#'     containing the values of measure points.
#'     \item \code{types}: The name(s) of the column(s) in \code{data_long}
#'     containing the type information of measure points (for different point
#'     color).
#'   }
#' @param estim_infos A list containing information about the estimations:
#'   \itemize{
#'     \item \code{col_mean}: The name of the column in \code{data_long}
#'     representing the mean value estimate.
#'     \item \code{col_ki_upp}: The name of the column in \code{data_long}
#'     representing the upper confidence band estimate.
#'     \item \code{col_ki_low}: The name of the column in \code{data_long}
#'     representing the lower confidence band estimate.
#'   }
#' @return A ggplot object.
#' @export
create_single_country_plot <- function(data_long,
                                       data_info,
                                       measure_info,
                                       estim_infos) {
  data_long2 <- data_long[data_long$country == data_info$name_country, ]
  measure_color <- "blue"
  ribbon_color <- "grey"
  
  out_plot <- ggplot2::ggplot(data_long2, ggplot2::aes(x = year)) +
    ggplot2::geom_point(ggplot2::aes_string(y = measure_info$vals,
                                            color = measure_info$types)) + 
    ggplot2::geom_line(ggplot2::aes_string(y = estim_infos$col_mean),
                                           color = measure_color) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymax = estim_infos$col_ki_upp,
                                             ymin = estim_infos$col_ki_low),
                         fill = ribbon_color, alpha = 0.25) +
    ggplot2::labs(title = data_info$name_country,
                  y = data_info$y_lab) +
    ggthemes::theme_tufte() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2")
  
  return(out_plot)
}
