% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/42_plot_measures.R
\name{create_me_plots_individual}
\alias{create_me_plots_individual}
\title{Generate Individual Measure Plots for Multiple Regressors}
\usage{
create_me_plots_individual(
  out_measures_info_KK,
  reg_names,
  settings = list(name_measure = "", plot_type = "base", plot_grid = c(4, 5), WITH_CI =
    TRUE)
)
}
\arguments{
\item{out_measures_info_KK}{A list where each element is a three-dimensional
array containing measure values for a specific regressor. Dimensions should
correspond to entities, time periods, and measure statistics (mean, upper,
and lower confidence intervals).}

\item{reg_names}{A vector of character strings specifying the names of
regressors for which plots will be generated. These names must correspond
to keys in `out_measures_info_KK`.}

\item{settings}{A list of settings for plot generation, which includes:
\itemize{
  \item{\code{name_measure}: }{The name of the measure being plotted, used
  for labeling.}
  \item{\code{plot_type}: }{Specifies the type of plot ("base" or "ggplot")
  .}
  \item{\code{plot_grid}: }{A vector specifying the layout of plots in
  terms of rows and columns (e.g., c(4, 5)).}
  \item{\code{WITH_CI}: }{Logical indicating whether to include confidence
  intervals in the plots. Defaults to `TRUE`.}
}}
}
\value{
The function implicitly returns a list of generated plot objects,
  although it is primarily used for its side effect of creating and saving
  plots.
}
\description{
Creates individual plots for each regressor within the specified measures
across different entities (e.g., countries) and time periods. This function
supports the inclusion of confidence intervals and allows customization of
plot settings.
}
