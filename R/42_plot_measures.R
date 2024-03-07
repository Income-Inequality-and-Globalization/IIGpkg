#' Generate Country-Specific Plots
#'
#' This function generates individual plots for each country based on Gini
#' coefficient data. It processes and reshapes the data, then creates and
#' arranges plots for each country.
#'
#' @param pth_data Path to the data file containing relevant data.
#' @param output_measure Array of estimated measure coefficient values as
#'   returned via [generate_measures_me()].
#' @param output_regs_grid Optional. A grid of regression outputs associated
#'   with measures, to be used for additional data points in the plots.
#'   Default is NULL, indicating no regression grid is used.
#' @param vars_to_use Vector of variable names to add as `geom_points` for
#'    comparison with a-posteriori estimation.
#' @param settings List of settings for data processing, including:
#' \itemize{
#' \item{\code{scale_transform}}{ Numeric value used to scale the measure values}
#' \item{\code{y_lab}}{ Label for the y-axis}
#' \item{\code{x_lab}}{ Label for the x-axis (optional)}
#' \item{\code{X_TRANSFORMED}}{ Boolean indicating whether to transform the x-axis
#'    values (optional)}
#' \item{\code{ADD_KI}}{ Boolean indicating whether to add confidence intervals
#'    to the plots (optional)}
#' }
#' @param name_measure Optional. Character string providing a name for the
#'   measure being plotted. This is used for naming within the plots if
#'   provided. Default is NULL.
#' @param PLOT Boolean indicating whether to actually plot the graphs or just
#'   return the data for plotting. Useful for non-interactive environments.
#'   Default is `TRUE`.
#'
#' @return A named list of two elements:
#'   \itemize{
#'     \item{\code{data_long}}{ Data in long format suitable for plotting}
#'     \item{\code{plot_objects}}{ A list of plot objects, one for each country,
#'        which can be rendered directly if PLOT is TRUE}
#'   }
#' @export
generate_country_plots <- function(pth_data,
                                   output_measure,
                                   output_regs_grid = NULL,
                                   vars_to_use,
                                   settings,
                                   name_measure = NULL,
                                   PLOT = TRUE) {
  names_measures <- get_names_measures(name_measure)
  data_to_plot <- get_data_measure_plots(output_measure,
                                         output_regs_grid,
                                         pth_data,
                                         vars_to_use,
                                         names_measures,
                                         settings)

  unique_countries <- unique(data_to_plot$country)
  list_plots_per_country <- list()
  # for (cntry in unique_countries[-c(length(unique_countries))]) {
  unique_countries <- unique_countries[-c(length(unique_countries))]
  for (cntry in unique_countries) {
    list_plots_per_country[[cntry]] <- create_single_country_plot(
      data_long = data_to_plot,
      name_country = cntry,
      plot_info = list(
        x_var = settings$x_var,
        y_lab = settings$y_lab,
        x_lab = settings$x_lab,
        X_TRN = settings$X_TRANSFORMED,
        ADD_KI = settings$ADD_KI
      ),
      measure_info = list(
        vals = "value",
        types = "type"
      ),
      estim_infos = names_measures
    )
  }
  if (PLOT) {
    layout_tkn <- get_layout_grid_from_sttgs(settings$grid_dim)
    gridExtra::grid.arrange(grobs = list_plots_per_country,
                            layout_matrix = layout_tkn,
                            as.table = FALSE)
  }
  return(list(data_long = data_to_plot,
              plot_objects = invisible(list_plots_per_country)))
}
get_layout_grid_from_sttgs <- function(sttgs_mfrow, transpose = TRUE) {
  out_grid <- matrix(seq_len(sttgs_mfrow[1] * sttgs_mfrow[2]),
                      nrow = sttgs_mfrow[1],
                      ncol = sttgs_mfrow[2],
                      byrow = transpose)
  return(out_grid)
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
                                 out_regs,
                                 pth_data,
                                 vars_to_use,
                                 names_measures,
                                 settings) {
  data_raw <- read.csv(file.path(pth_data)) %>% tibble::as_tibble()
  
  scale_values <- get_default_scale_measre(settings$scale_measure)
  reg_name_tkn <- settings$x_var

  measure_mean_01 <- get_measure_vals(out_measures, "mean", scale_values)
  measure_ki_l_02 <- get_measure_vals(out_measures, "ki_low", scale_values)
  measure_ki_u_03 <- get_measure_vals(out_measures, "ki_upp", scale_values)
  data_raw[names_measures$mean]   <- measure_mean_01
  data_raw[names_measures$ki_low] <- measure_ki_l_02
  data_raw[names_measures$ki_upp] <- measure_ki_u_03

  country_list <- unique(data_raw$country)
  year_list <- unique(data_raw$year)
  data_regs <- get_tibble_regs_grid(out_regs,
                                    reg_name_tkn,
                                    country_list,
                                    year_list)

  vars_required <- c(unlist(names_measures, use.names = FALSE),
                     vars_to_use,
                     reg_name_tkn)
  data_raw <- data_raw %>%
    dplyr::select("country", "year", tidyr::all_of(vars_required))
  if (!is.null(data_regs)) {
    data_raw <- data_raw %>%
      get_sorted_reg("country", reg_name_tkn) %>%
      dplyr::full_join(data_regs, by = c("country", "year"))
    data_raw$year <- NULL
  }
  return(data_raw)
}
get_measure_vals <- function(out_msr, name_var, scl_val = 1) {
  as.vector(t(out_msr[, , name_var])) * scl_val
}
get_default_scale_measre <- function(scl_msr) {
  if (is.null(scl_msr)) return(1)
  scl_msr
}
get_sorted_reg <- function(df, name_country_col, reg_to_sort) {
  unq_countries   <- unique(df[[name_country_col]])
  reg_sorted <- numeric(nrow(df))  # Initialize a vector for sorted values

  for (cntry in unq_countries) {
    id_rows <- which(df[[name_country_col]] == cntry)
    reg_sorted[id_rows] <- sort(df[[reg_to_sort]][id_rows], na.last = TRUE)
  }
  df[[reg_to_sort]] <- reg_sorted
  return(df)
}

get_tibble_regs_grid <- function(regs_grid,
                                 reg_name = NULL,
                                 country_list,
                                 year_list) {
  if (is.null(regs_grid) || is.null(reg_name)) return(NULL)
  reg_names_avail <- dimnames(regs_grid)[[3]]
  stopifnot(reg_name %in% reg_names_avail)

  id_regs  <- which(grepl(reg_name, dimnames(regs_grid)[[1]]))
  regs_out <- tibble::as_tibble(t(regs_grid[id_regs, , reg_name]))
  colnames(regs_out) <- country_list
  regs_out$year <- year_list
  reg_name_trns <- paste0(reg_name, "_transformed")
  regs_out <- regs_out %>%
    tidyr::pivot_longer(cols = !"year",
                        names_to = "country",
                        values_to = reg_name_trns)
  regs_out <- regs_out %>%
    dplyr::select("country", "year", reg_name_trns) %>%
    dplyr::arrange(.data$country, .data$year)

  return(regs_out)
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
#' @param data_regs_grid an array of dimension `(NN x num_par) x TT x KK` where
#'   the third  dimension gives the number of regressors (the second time and
#'   the first a is cross section )
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
                                   output_regs_grid,
                                   pth_data_gmm,
                                   vars_to_use,
                                   names_to_use,
                                   settings) {
  data_processed <- data_plot_processing(data_estim,
                                         output_regs_grid,
                                         pth_data_gmm,
                                         vars_to_use,
                                         names_to_use,
                                         settings)
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
#' @param name_country character giving the name of the country for the plot
#' @param plot_info A list containing information about the data and plot:
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
                                       name_country,
                                       plot_info,
                                       measure_info,
                                       estim_infos) {
  data_long2 <- data_long[data_long$country == name_country, ]
  measure_color <- "blue"
  ribbon_color <- "darkgrey"

  x_reg <- plot_info$x_var
  if (is.null(x_reg)) x_reg <- "year"
  if (!is.null(x_reg)) {
    if (isTRUE(plot_info$X_TRN)) {
      paste0(x_reg, "_transformed")
    }
  }
  out_plot <- ggplot2::ggplot(data_long2, ggplot2::aes_string(x = x_reg)) +
    ggplot2::geom_point(ggplot2::aes_string(y = measure_info$vals,
                                            color = measure_info$types)) +
    ggplot2::geom_line(ggplot2::aes_string(y = estim_infos$mean),
                                           color = measure_color)
  if (isTRUE(plot_info$ADD_KI)) {
    out_plot <- out_plot +
      ggplot2::geom_ribbon(ggplot2::aes_string(ymax = estim_infos$ki_upp,
                                               ymin = estim_infos$ki_low),
                           fill = ribbon_color,
                           alpha = 0.5)
  }
  out_plot <- out_plot +
    ggplot2::labs(title = name_country,
                  y = plot_info$y_lab) +
    ggthemes::theme_tufte() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2")

  return(out_plot)
}
#' Generate Individual Measure Plots for Multiple Regressors
#'
#' This function creates individual plots for each regressor within the
#' specified measures. It iterates over each measure, generating plots for each
#' point in the measure's dimension (e.g., time series data for multiple
#' countries). The function supports the inclusion of confidence intervals
#' based on the settings provided.
#'
#' @param out_measures_info_KK A list where each element is a three-dimensional
#'   array of measure values for a specific regressor. The dimensions should
#'   correspond to different entities (e.g., countries), time periods, and
#'   measure statistics (mean, upper, and lower confidence intervals).
#' @param reg_names Vector of character strings specifying the names of the
#'   regressors for which plots will be generated. These names must correspond
#'   to keys in `out_measures_info_KK`.
#' @param settings List of settings for plot generation, including:
#'   \itemize{
#'     \item `mfrow`: A vector specifying the layout of plots in terms of rows
#'     and columns.
#'     \item `WITH_CI`: Logical indicating whether to include confidence
#'     intervals in the plots. Defaults to TRUE.
#'   }
#'
#' @return Generates plots for each specified regressor and measure combination
#'   but does not explicitly return any value.
#' @export
create_me_plots_individual <- function(
    out_measures_info_KK,
    reg_names,
    settings =
        list(name_measure = "",
             plot_grid = c(4, 5),
             WITH_CI = TRUE)) {
  num_regs_me  <- length(reg_names)
  info_on_plot <- dimnames(out_measures_info_KK[[reg_names[1]]])

  NN <- length(info_on_plot[[1]])
  TT <- length(info_on_plot[[2]])

  y_lab     <- settings$name_measure
  WITH_CI   <- settings$WITH_CI
  plot_grid <- settings$plot_grid

  pdf_file_name <- paste0("03_",settings$name_measure, "_CI_", ".pdf")
  pdf(pdf_file_name, width = 11, height = 8.5)
  par(mfrow = plot_grid)
  for (nn in seq_len(NN)) {
    for (kk in seq_len(num_regs_me)) {

      par(mfrow = plot_grid)
      for (tt in 2:TT) {
        vals_to_plot <- out_measures_info_KK[[reg_names[kk]]][nn, tt, , ]
        min_max <- get_min_max_y_scale(vals_to_plot)
        get_single_plot_me(
          vals_to_plot,
          settings = list(WITH_CI = WITH_CI,
                          title = paste0("year: ", info_on_plot[[2]][tt]),
                          type = "plot",
                          ggplot = TRUE,
                          line_col = "black",
                          y_lab = y_lab,
                          min_max = min_max)
        )
      }
      mtext(paste0("Country: ",
                   info_on_plot[[1]][nn],
                   " -- reg: ",
                   reg_names[kk]),
            side = 3,
            line = -1.5,
            outer = TRUE)
      # Close PDF device to save the plot
    }
  }
  dev.off()
  par(mfrow = c(1, 1))
}
#' Generate Time Series Plots for Multiple Regressors
#'
#' This function creates time series plots for each regressor across all
#' entities within the specified measures. Unlike `create_me_plots_individual`,
#' which generates plots for each time period separately, this function
#' integrates time periods into a single plot per entity, facilitating a
#' continuous time series visualization. Confidence intervals are not included
#' in this visualization, focusing instead on the trend over time.
#'
#' @param out_measures_info_KK A list where each element is a three-dimensional
#'   array of measure values for a specific regressor. The dimensions should
#'   correspond to different entities (e.g., countries), time periods, and
#'   measure statistics (mean values only, as confidence intervals are not
#'   displayed).
#' @param reg_names Vector of character strings specifying the names of the
#'   regressors for which plots will be generated. These names must correspond
#'   to keys in `out_measures_info_KK`.
#' @param settings List of settings for plot generation, including:
#'   \itemize{
#'     \item `mfrow`: A vector specifying the layout of plots in terms of rows
#'     and columns. This controls the arrangement of time series plots for
#'     multiple entities.
#'     \item `WITH_CI`: Logical indicating whether to include confidence
#'     intervals in the plots. In this function, it is effectively ignored as
#'     confidence intervals are not supported.
#'   }
#'
#' @return Generates time series plots for each specified regressor and entity
#'   combination but does not explicitly return any value. The function
#'   emphasizes trend visualization over time for each entity across all
#'   specified regressors.
#' @export

create_me_plots_time_series <- function(
    out_measures_info_KK,
    reg_names,
    settings =
        list(name_measure = "",
             plot_type = "base",
             plot_grid = c(3, 5))) {
  num_regs_me  <- length(reg_names)
  info_on_plot <- dimnames(out_measures_info_KK[[reg_names[1]]])

  NN <- length(info_on_plot[[1]])
  TT <- length(info_on_plot[[2]])
  title_nn <- substr(info_on_plot[[1]], 1, 5)
  plot_type <- settings$plot_type

  col_off <- 5
  color_palette <- scales::col_numeric(palette = "Blues",
                                       domain = c(1, TT + col_off))
  y_lab     <- settings$name_measure
  plot_grid <- settings$plot_grid
  out_plot  <- NULL
  out_plot_list <- list()
  iter_plot <- 1

  if (plot_type == "base") {
    par(mfcol = settings$plot_grid)
    ggplot <- FALSE
  } else if (plot_type == "ggplot") {
    ggplot <- TRUE
  }
  for (nn in seq_len(NN)) {
    for (kk in seq_len(num_regs_me)) {
      base_ggplot <- NULL
      vals_to_plot <- out_measures_info_KK[[reg_names[kk]]][nn, , , ]
      min_max <- get_min_max_y_scale(vals_to_plot[, , 1])
      base_ggplot <- get_single_plot_me(vals_to_plot[1, , ],
                                        settings = list(
                                          WITH_CI = FALSE,
                                          type = "plot",
                                          ggplot = ggplot,
                                          y_lab = y_lab,
                                          x_lab = reg_names[kk],
                                          title = title_nn[nn],
                                          min_max = min_max,
                                          line_col = color_palette(1 + col_off))
                                        )
      for (tt in 2:TT) {
        # Compute color based on time point
        line_color_tkn <- color_palette(tt + col_off)
        base_ggplot <- get_single_plot_me(vals_to_plot[tt, , ],
                                          settings = list(
                                            WITH_CI = FALSE,
                                            type = "line",
                                            ggplot = ggplot,
                                            base_ggplot = base_ggplot,
                                            line_col = line_color_tkn)
                                          )
      }
      out_plot_list[[iter_plot]] <- base_ggplot
      iter_plot <- iter_plot + 1
    }
    if (plot_type == "base" && nn %% plot_grid[2] == 0) {
      mtext(paste0("Marginal effect on: ",
                   settings$name_measure),
            side = 3,
            line = -1.5,
            outer = TRUE)
    }
  }
  if (plot_type == "base") {
    par(mfrow = c(1, 1))
  }
  if (plot_type == "ggplot") {
    layout_tkn <- get_layout_grid_from_sttgs(plot_grid, transpose = TRUE)
    out_plot <- gridExtra::marrangeGrob(
      grobs = out_plot_list,
      layout_matrix = layout_tkn,
      as.table = FALSE)
  }
  return(invisible(list(out_plot_list = out_plot_list, out_plot = out_plot)))
}
#' Plot Single Measure with Optional Confidence Intervals, Custom Title, and
#' Type
#'
#' This function generates a plot for a single measure, with the option to
#' include confidence intervals for a detailed analysis, and allows
#' customization of the plot type. It visualizes the measure across a specified
#' dimension, such as time. The function supports customization including plot
#' titles, the decision to include confidence intervals, and the type of plot
#' (line or plot).
#'
#' @param vals_to_plot Matrix where the first column contains the measure values
#'   to be plotted. If `WITH_CI` is TRUE within the settings list, the second
#'   and third columns must contain the upper and lower confidence intervals,
#'   respectively.
#' @param settings A list of settings for the plot, including:
#'   \itemize{
#'     \item `WITH_CI`: Logical indicating whether to include confidence
#'     intervals. If FALSE, only the measure values are plotted. If TRUE, dashed
#'     lines represent the confidence intervals. Cannot be TRUE when type is
#'     "line".
#'     \item `title`: The title of the plot. A string that will be displayed as
#'     the main title of the plot.
#'     \item `type`: Character indicating the type of plot. Accepts "plot" for a
#'     basic plot or "line" for a line plot. The option "line" cannot have
#'     confidence intervals.
#'   }
#'
#' @return A plot is generated with the specified parameters but not explicitly
#'   returned. The plot includes measure values and, optionally, confidence
#'   intervals, along with a custom title and specified plot type.
#'
#' @export
get_single_plot_me <- function(
  vals_to_plot,
  settings =
    list(WITH_CI = FALSE,
         y_lab = "",
         x_lab = "",
         title = "",
         line_col = "black",
         min_max = NULL,
         type = "",
         ggplot = FALSE,
         base_ggplot = NULL)
) {
  WITH_CI <- settings$WITH_CI
  title   <- settings$title
  type    <- settings$type
  ggplot  <- settings$ggplot
  base_p  <- settings$base_ggplot
  y_lab   <- settings$y_lab
  x_lab   <- settings$x_lab
  min_max <- settings$min_max
  line_col <- settings$line_col
  stopifnot(`Argument 'type' must unknown.` = type %in% c("plot", "line"))
  stopifnot(`Arg. 'ggplot' is not logical.` = is.logical(ggplot))
  if (ggplot) {
    plot_data <- data.frame(x = seq_along(vals_to_plot[, 1]),
                            y = vals_to_plot[, 1],
                            upper = vals_to_plot[, 2],
                            lower = vals_to_plot[, 3])
    if (type == "plot") {
      me_plot <- plot_data %>%
        ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_line(color = line_col) +
        # ggplot2::scale_color_manual(values = "#DEEBF7") +
        ggplot2::labs(title = settings$title,
                      y = settings$y_lab,
                      x = settings$x_lab) +
        ggplot2::theme_minimal()
    } else if (type == "line") {
      stopifnot("Arg. 'base_plot' must be a ggplot object." = !is.null(base_p))
      me_plot <- base_p + ggplot2::geom_line(ggplot2::aes(y =  plot_data$y),
                                             color = line_col)
    }
    # Conditionally add CIs with ribbon if WITH_CI is TRUE
    if (settings$WITH_CI) {
      if (type == "line") stop("Cannot have CI bands with type = 'line'.")
      me_plot <- me_plot +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper),
                                          fill = "grey80",
                                          col = "darkred",
                                          alpha = 0.5)
    }
    return(me_plot)
  } else {
    if (isTRUE(WITH_CI)) {
      if (type == "line") stop("Cannot have CI bands with type = 'line'.")
      mean_to_plot <- vals_to_plot[, 1]
      ki_upp       <- vals_to_plot[, 2]
      ki_low       <- vals_to_plot[, 3]

      plot(mean_to_plot, type = "l",
           ylim = c(min_max[1], min_max[2]),
           main = title,
           ylab = y_lab,
           xlab = x_lab,
           col = line_col)
      lines(ki_upp, col = "darkred", lty = "dashed")
      lines(ki_low, col = "darkred", lty = "dashed")
    } else {
      if (type == "plot") {
        plot(vals_to_plot[, 1],
             type = "l",
             ylim = c(min_max[1], min_max[2]),
             main = title,
             ylab = y_lab,
             xlab = x_lab,
             col = line_col)
      } else if (type == "line") {
        lines(vals_to_plot[, 1], col = line_col)
      }
    }
  me_plot <- NULL
  return(invisible(me_plot))
  }
}
get_min_max_y_scale <- function(x) {
  min_y <- min(x)
  max_y <- max(x)
  return(c(min_y, max_y))
}