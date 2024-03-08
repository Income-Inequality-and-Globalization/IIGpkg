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
#' \item{\code{scale_transform}}{ Numeric value used to scale the measure
#' values}
#' \item{\code{y_lab}}{ Label for the y-axis}
#' \item{\code{x_lab}}{ Label for the x-axis (optional)}
#' \item{\code{X_TRANSFORMED}}{ Boolean indicating whether to transform the
#' x-axis values (optional)}
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
#'     \item{\code{plot_objects}}{ A list of plot objects which can be rendered
#'     directly if PLOT is TRUE}
#'   }
#' @export
generate_country_plots <- function(pth_data,
                                   output_measure,
                                   output_regs_grid = NULL,
                                   vars_to_use,
                                   settings,
                                   name_measure = NULL,
                                   PLOT = TRUE) {
  nm_pdf   <- "01_fitted_mes.pdf"
  grid_dim <- settings$grid_dim
  out_TT <- output_measure$info_TT
  out_KK <- output_measure$info_KK
  names_measures <- get_names_measures(name_measure)
  data_to_plot <- get_data_measure_plots(out_TT,
                                         output_regs_grid,
                                         pth_data,
                                         vars_to_use,
                                         names_measures,
                                         settings)

  unique_countries <- unique(data_to_plot$country)
  list_plots_per_country <- list()
  # unique_countries <- unique_countries[-c(length(unique_countries))]
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
  list_plots_me <- create_me_plots_time_series(
    out_KK,
    output_regs_grid,
    regs_to_use,
    RENDER_PLOT = FALSE,
    NO_TITLE = TRUE,
    settings = list(name_measure = name_measure,
                    plot_type = "ggplot"))[["out_plot_list"]]
  list_plots <- get_sorted_plot_list_all(list_plots_per_country, list_plots_me)
  if (PLOT) {
    grid_dim_cnt <- c(1, grid_dim[2])
    grid_dim_mes <- c(grid_dim[1] - 1, grid_dim[2])
    lmat  <- rbind(
      get_layout_grid_from_sttgs(grid_dim_cnt),
      get_layout_grid_from_sttgs(grid_dim_mes, transpose = FALSE) + grid_dim_cnt[2])
    glist <- lapply(list_plots, ggplot2::ggplotGrob)
    ggplot2::ggsave(
      nm_pdf,
      gridExtra::marrangeGrob(glist,
                              layout_matrix = lmat,
                              nrow = grid_dim[1],
                              ncol = grid_dim[2]),
      width = 29.7, height = 21, units = "cm")
  }
  return(list(data_long = data_to_plot,
              plot_objects = list_plots))
}
get_layout_grid_from_sttgs <- function(sttgs_mfrow, transpose = TRUE) {
  out_grid <- matrix(seq_len(sttgs_mfrow[1] * sttgs_mfrow[2]),
                      nrow = sttgs_mfrow[1],
                      ncol = sttgs_mfrow[2],
                      byrow = transpose)
  return(out_grid)
}
get_sorted_plot_list_all <- function(plot_list_cnt, plot_list_mes) {
  num_plot_cnt <- length(plot_list_cnt)
  num_plot_mes <- length(plot_list_mes)
  num_half_cnt <- num_plot_cnt / 2
  num_half_mes <- num_plot_mes / 2
  c(head(plot_list_cnt, num_half_cnt),
    head(plot_list_mes, num_half_mes),
    tail(plot_list_cnt, num_half_cnt),
    tail(plot_list_mes, num_half_mes))
}
#' Read and Preprocess Data for Plotting
#'
#' This function reads data from a specified path, preprocesses it by scaling
#' measure-coefficient values, and optionally merges it with regression data.
#' The processed data is then ready for plotting or further analysis. This
#' includes scaling the measure values, merging additional regression data, and
#' selecting relevant variables for visualization.
#'
#' @param out_measures A three-dimensional array of estimated
#'   measure-coefficient values. The dimensions are expected to correspond to
#'   entities (e.g., countries), time periods, and measure statistics (mean,
#'   lower confidence band, upper confidence band). This data is used to append
#'   new variables to the dataset.
#' @param out_regs An optional data frame or matrix containing regression data
#'   to be merged with the measure data. Each row should correspond to an
#'   observation, and columns should include the regression variables of
#'   interest. If NULL, no regression data is merged.
#' @param pth_data Path to the CSV data file containing the raw data set. The
#'   file is read into R and processed according to the specified parameters.
#' @param vars_to_use A vector of variable names from the raw data set that
#'   should be retained for analysis and plotting. This allows for the inclusion
#'   of specific variables of interest beyond the generated measure variables.
#' @param names_measures A named list or vector specifying the column names to
#'   be created in the dataset for each type of measure (mean, lower confidence
#'   band, upper confidence band). These names are dynamically used to assign
#'   measure values to the appropriate columns in the processed dataset.
#' @param settings A list of settings for data processing, including:
#'   \itemize{
#'     \item{\code{scale_measure}:}{ A numeric value used to scale the measure
#'     values. Defaults to 1 if not specified, indicating no scaling. This is
#'     applied uniformly to mean, lower, and upper confidence band values.}
#'     \item{\code{x_var}:}{ The name of the variable used for regression data
#'     merging and/or plotting. This should correspond to a variable in both the
#'     raw data and the regression data (if provided).}
#'     \item Other settings as required for processing, such as transformation
#'     flags or additional variable names to include.
#'   }
#'
#' @return A tibble containing the processed data. This includes the original
#'   variables specified by `vars_to_use`, new measure variables as named by
#'   `names_measures`, and optionally merged regression data. Variables not
#'   specified for retention or addition are omitted from the returned tibble.
data_plot_processing <- function(out_measures,
                                 out_regs,
                                 pth_data,
                                 vars_to_use,
                                 names_measures,
                                 settings) {
  data_raw <- read.csv(file.path(pth_data)) %>% tibble::as_tibble()
  data_raw$est_mean <- data_raw$est_mean * 100

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
#' This function generates a ggplot object for a specified country, plotting
#' measure values across time. It visualizes the data in a way that highlights
#' changes in measure values over the years, using points and lines, and
#' optionally includes a confidence interval ribbon if specified.
#'
#' @param data_long The data in long format, containing variables such as
#'   country, year, and the measures to be plotted.
#' @param name_country A character string giving the name of the country for
#'   which the plot will be generated. This specifies the subset of data to use
#'   for the plot.
#' @param plot_info A list containing settings for the plot, including:
#'   \itemize{
#'     \item{x_var}{The variable name from `data_long` to be used as the x-axis.
#'     Typically, this would be the year or time variable.}
#'     \item{y_lab}{Label for the y-axis, typically describing the measure being
#'     plotted.}
#'     \item{x_lab}{Label for the x-axis, typically "Year". Optional.}
#'     \item{X_TRN}{Logical indicating whether the x-axis variable should be
#'     transformed. Optional.}
#'     \item{ADD_KI}{Logical indicating whether to add a confidence interval
#'     ribbon to the plot.}
#'   }
#' @param measure_info A list specifying the columns in `data_long` related to
#'   the measures:
#'   \itemize{
#'     \item{vals}{The name of the column containing the values of the measure
#'     points to be plotted.}
#'     \item{types}{The name of the column indicating the type of measure point
#'     (e.g., for differentiating points by color). Optional.}
#'   }
#' @param estim_infos A list detailing columns in `data_long` for estimation
#'   information:
#'   \itemize{
#'     \item{mean}{The column name representing the mean value estimate of the
#'     measure.}
#'     \item{ki_upp}{The column name for the upper confidence band estimate.
#'     Required if `ADD_KI` is TRUE.}
#'     \item{ki_low}{The column name for the lower confidence band estimate.
#'     Required if `ADD_KI` is TRUE.}
#'   }
#' @return A `ggplot` object representing the specified country's measure values
#'   over time, with optional confidence intervals.
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
#' Creates individual plots for each regressor within the specified measures
#' across different entities (e.g., countries) and time periods. This function
#' supports the inclusion of confidence intervals and allows customization of
#' plot settings.
#'
#' @param out_measures_info_KK A list where each element is a three-dimensional
#'   array containing measure values for a specific regressor. Dimensions should
#'   correspond to entities, time periods, and measure statistics (mean, upper,
#'   and lower confidence intervals).
#' @param reg_names A vector of character strings specifying the names of
#'   regressors for which plots will be generated. These names must correspond
#'   to keys in `out_measures_info_KK`.
#' @param settings A list of settings for plot generation, which includes:
#'   \itemize{
#'     \item{\code{name_measure}: }{The name of the measure being plotted, used
#'     for labeling.}
#'     \item{\code{plot_type}: }{Specifies the type of plot ("base" or "ggplot")
#'     .}
#'     \item{\code{plot_grid}: }{A vector specifying the layout of plots in
#'     terms of rows and columns (e.g., c(4, 5)).}
#'     \item{\code{WITH_CI}: }{Logical indicating whether to include confidence
#'     intervals in the plots. Defaults to `TRUE`.}
#'   }
#'
#' @return The function implicitly returns a list of generated plot objects,
#'   although it is primarily used for its side effect of creating and saving
#'   plots.
#' @export

create_me_plots_individual <- function(
    out_measures_info_KK,
    reg_names,
    settings =
        list(name_measure = "",
             plot_type = "base",
             plot_grid = c(4, 5),
             WITH_CI = TRUE)) {
  num_regs_me  <- length(reg_names)
  info_on_plot <- dimnames(out_measures_info_KK[[reg_names[1]]])

  NN <- length(info_on_plot[[1]])
  TT <- length(info_on_plot[[2]])

  y_lab     <- settings$name_measure
  WITH_CI   <- settings$WITH_CI
  plot_type <- settings$plot_type
  plot_grid <- settings$plot_grid
  out_plot  <- NULL
  out_plot_list <- list()
  iter_plot <- 1

  pdf_file_name <- paste0("03_",settings$name_measure, "_CI_", ".pdf")
  pdf(pdf_file_name, width = 11, height = 8.5)
  if (plot_type == "base") {
    par(mfcol = plot_grid)
    ggplot <- FALSE
  } else if (plot_type == "ggplot") {
    ggplot <- TRUE
  }
  for (nn in seq_len(NN)) {
    for (kk in seq_len(num_regs_me)) {
      for (tt in 2:TT) {
        vals_to_plot <- out_measures_info_KK[[reg_names[kk]]][nn, tt, , ]
        min_max <- get_min_max_y_scale(vals_to_plot)
        out_plot_list[[iter_plot]] <- get_single_plot_me(
          vals_to_plot,
          settings = list(WITH_CI = WITH_CI,
                          title = paste0("year: ", info_on_plot[[2]][tt]),
                          type = "plot",
                          ggplot = ggplot,
                          line_col = "black",
                          y_lab = y_lab,
                          min_max = min_max)
        )
        iter_plot <- iter_plot + 1
      }
      if (plot_type == "base") {
        mtext(paste0("Country: ",
                     info_on_plot[[1]][nn],
                     " -- reg: ",
                     reg_names[kk]),
              side = 3,
              line = -1.5,
              outer = TRUE)
      }
    }
  }
  # Close PDF device to save the plot
  dev.off()
  par(mfrow = c(1, 1))
  return(invisible(list(out_plot_list = out_plot_list, out_plot = out_plot)))
}
#' Generate Time Series Plots for Multiple Regressors
#'
#' Creates time series plots for each regressor across all entities within the
#' specified measures, facilitating a continuous time series visualization.
#' Unlike `create_me_plots_individual`, this function integrates time periods
#' into a single plot per entity, and confidence intervals are not included.
#'
#' @param out_measures_info_KK A list where each element is a three-dimensional
#'   array of measure values for a specific regressor. The dimensions should
#'   correspond to different entities (e.g., countries), time periods, and
#'   measure statistics (mean values only).
#' @param reg_names A vector of character strings specifying the names of the
#'   regressors for which plots will be generated. These names must correspond
#'   to keys in `out_measures_info_KK`.
#' @param settings A list of settings for plot generation, including:
#'   \itemize{
#'     \item{\code{name_measure}: }{The measure's name, used for labeling.}
#'     \item{\code{plot_type}: }{Specifies the plot's type ("base" or "ggplot")
#'     .}
#'     \item{\code{plot_grid}: }{A vector specifying the layout of plots in
#'     terms of rows and columns.}
#'     \item{\code{WITH_CI}: }{This parameter is ignored in this function as
#'     confidence intervals are not included in the visualization.}
#'   }
#'
#' @return Implicitly returns a list of plot objects or a composite graphic,
#'   depending on the `plot_type`. Primarily used for its side effect of
#'   plotting.
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
#' Plot Single Measure with Confidence Intervals, Custom Title, and Type
#'
#' Generates a plot for a single measure, allowing for confidence intervals,
#' customization of plot type, and visualizing the measure across a specified
#' dimension. Supports both base R plots and ggplot2, with options for plot
#' customization.
#'
#' @param vals_to_plot A matrix where columns represent the measure values to be
#'   plotted, and, if `WITH_CI` is TRUE, the second and third columns must
#'   contain the upper and lower confidence intervals, respectively.
#' @param settings A list of settings for the plot, which includes:
#'   \itemize{
#'     \item{\code{WITH_CI} :}{Logical indicating whether to include confidence
#'     intervals. If `FALSE`, only measure values are plotted. If `TRUE`, and
#'     type is "plot", confidence intervals are shown as dashed lines.}
#'     \item{\code{title} :}{The title of the plot.}
#'     \item{\code{type} :}{The type of plot ("plot" for a basic plot or "line"
#'     for a line plot). "line" type cannot include confidence intervals.}
#'     \item{\code{ggplot} :}{Logical indicating if ggplot2 should be used.
#'     Defaults to FALSE.}
#'     \item{\code{y_lab} :}{Y-axis label.}
#'     \item{\code{x_lab} :}{X-axis label.}
#'     \item{\code{line_col} :}{Color of the line or points.}
#'     \item{\code{min_max} :}{Vector specifying the minimum and maximum y-axis
#'     values.}
#'     \item{\code{base_ggplot} :}{An existing ggplot object to add the line or
#'     points to, required if `type` is "line" and `ggplot` is TRUE.}
#'   }
#'
#' @return Generates and potentially returns a plot object based on the
#'   specified parameters, including measure values and optional confidence
#'   intervals.
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