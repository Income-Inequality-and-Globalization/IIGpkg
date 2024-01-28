#' Generate Country-Specific Plots
#'
#' This function generates individual plots for each country based on Gini
#' coefficient data. It processes and reshapes the data, then creates and
#' arranges plots for each country.
#'
#' @param pth_data path to the data file containing relevant data.
#' @param output_measure Array of estimated measure coefficient values as
#'   returned via [generate_measures_me()].
#' @param names_measures character giving the name of the measure; can be left
#'   `NULL` in which case the internal data set is getting a default placeholder
#'   for the variable names
#' @param vars_to_use Vector of variable names to add as `geom_points` for
#'    comparison with a-posteriori estimation.
#' @param settings List of settings for data processing, including:alpha
#' \itemize{
#' \item{scale_transform}{Numeric value used to scale the measure values}
#' \item{y_lab}{Label for the y-axis}
#' }
#'
#' @return A named list of three elements:
#'    \itemize{
#'       \item{`data_long`}{long format data}
#'       \item{`plot_objects`}{grid of plots, one for each country}
#'      }
#' @export
generate_country_plots <- function(pth_data,
                                   output_measure,
                                   vars_to_use,
                                   settings,
                                   name_measure = NULL,
                                   PLOT = TRUE) {
  name_measure <- paste0(name_measure, "_")
  names_measures <- list(
    mean = paste0(name_measure, "mean"),
    ki_low = paste0(name_measure, "ki_low"),
    ki_upp = paste0(name_measure, "ki_upp")
  )
  data_to_plot <- get_data_measure_plots(output_measure,
                                         pth_data,
                                         vars_to_use,
                                         names_measures,
                                         settings)

  unique_countries <- unique(data_to_plot$country)
  list_plots_per_country <- list()
  for (country in unique_countries[-c(length(unique_countries))]) {
    list_plots_per_country[[country]] <- create_single_country_plot(
      data_long = data_to_plot,
      data_info = list(
        name_country = country,
        y_lab = settings$y_lab
      ),
      measure_info = list(
        vals = "value",
        types = "type"
      ),
      estim_infos = names_measures
    )
  }
  if (PLOT) {
    gridExtra::grid.arrange(grobs = list_plots_per_country, as.table = FALSE)
  }
  return(list(data_long = data_to_plot,
              plot_objects = invisible(list_plots_per_country)))
}
#' Read and Preprocess Data
#'
#' Reads data from a specified `path` and `filename`, and preprocesses it by
#' removing specified variables, and adding output measures from a suitable
#' array (`out_measures`).
#'
#' @param out_measures an array of estimated measure-coefficient values:
#'   rownmaes give the country (in the same order as the dataset), colnames the
#'   time, and the third dimension has the Bayesian a-posteriori mean, lower and
#'   upper confidence band values. The measure values can be Ginis, or the
#'   expectaion as estimated from the Bayesian parameter samples a-posteriori.
#' @param pth_data The path to the data file which appended by content of the
#'   of the first argument (ie. transform the arrary to the data format and
#'   add its content as new variables)
#' @param names_measures A named list or vector specifying the column names to
#'   be created in the dataset for each type of measure. It must include names
#'   for "mean", "ki_low", and "ki_upp", corresponding to the Bayesian
#'   a-posteriori mean, lower confidence band, and upper confidence band values,
#'   respectively. The names are used to dynamically create and assign these
#'   measures to the appropriate columns in the dataset. For example,
#'   names_measures could be set as `list(mean = "gini_mean", ki_low =
#'   "gini_ki_low", ki_upp = "gini_ki_upp")`.
#' @param scale_values A numeric value used to scale the measure values
#'   extracted from `out_measures`. This scaling factor is applied to the mean,
#'   lower, and upper confidence band values. It allows for adjusting the
#'   measure values according to the desired scale (e.g., converting proportions
#'   to percentages). Default is 1 (no scaling).
#'
#' @return A tibble with the specified variable removed, and the output measure
#'   variables added
data_plot_processing <- function(out_measures,
                                 pth_data,
                                 names_measures,
                                 scale_values = 1) {
  data_raw <- read.csv(file.path(pth_data)) %>% tibble::as_tibble()

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
#' variables to be reshaped, making it versatile for various types of data
#' beyond just Gini coefficients.
#'
#' @param data_raw The raw data as a tibble.
#' @param var_to_use A character vector of variable names to be used in the
#'   transformation. These variables will be reshaped from wide to long format.
#' @param nms_to A string specifying the name of the new column in the long
#'   format that will contain the variable names from `var_to_use`.
#' @param vals_to A string specifying the name of the new column in the long
#'   format that will contain the values from the variables in `var_to_use`.
#'
#' @return A tibble in long format with columns as specified by `nms_to` and
#'   `vals_to`.
get_data_plot_long <- function(data_raw, var_to_use, nms_to, vals_to) {
  data_long <- data_raw %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(var_to_use),
                        names_to = nms_to,
                        values_to = vals_to)
  return(data_long)
}

#' Generate Plots for Data Measures
#'
#' This function serves as a wrapper for `data_plot_processing` and
#' `get_data_plot_long`, facilitating the generation of measure-based plots. It
#' first processes the data by incorporating estimated measures and then
#' reshapes it into a long format suitable for plotting.
#'
#' @param data_estim An array of estimated measure-coefficient values,
#' structured similarly to the `out_measures` parameter in
#'   `data_plot_processing`. This array is used to add new variables to the
#'   dataset based on the estimated measures.
#' @param pth_data_gmm The path to the data file, which is used by the
#'   `data_plot_processing` function to read and preprocess the data.
#' @param vars_to_use A character vector specifying the variable names to be
#'   included in the long format data. These are the variables that will be
#'   reshaped and used for plotting.
#' @param names_to_use A named list or vector specifying the column names to be
#'   created in the dataset for each type of measure, used in the
#'   `data_plot_processing` function.
#' @param settings A list containing settings for data processing. This should
#'   include at least `scale_transform`, a numeric value used to scale the
#'   measure values in `data_plot_processing`.
#'
#' @return A tibble in long format, suitable for plotting, containing the
#'   reshaped data with added measure variables.
#'
#' @export
get_data_measure_plots <- function(data_estim,
                                   pth_data_gmm,
                                   vars_to_use,
                                   names_to_use,
                                   settings) {
  data_processed <- data_plot_processing(data_estim,
                                         pth_data_gmm,
                                         names_to_use,
                                         settings$scale_transform)
  data_long <- get_data_plot_long(data_processed,
                                  vars_to_use,
                                  "type",
                                  "value")
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
#'     \item \code{mean}: The name of the column in \code{data_long}
#'     representing the mean value estimate.
#'     \item \code{ki_upp}: The name of the column in \code{data_long}
#'     representing the upper confidence band estimate.
#'     \item \code{ki_low}: The name of the column in \code{data_long}
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
    ggplot2::geom_line(ggplot2::aes_string(y = estim_infos$mean),
                                           color = measure_color) +
    ggplot2::geom_ribbon(ggplot2::aes_string(ymax = estim_infos$ki_upp,
                                             ymin = estim_infos$ki_low),
                         fill = ribbon_color, alpha = 0.25) +
    ggplot2::labs(title = data_info$name_country,
                  y = data_info$y_lab) +
    ggthemes::theme_tufte() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2")

  return(out_plot)
}
