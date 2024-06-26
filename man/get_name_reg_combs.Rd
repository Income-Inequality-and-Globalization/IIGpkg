% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_helper_data_cleaning.R
\name{get_name_reg_combs}
\alias{get_name_reg_combs}
\title{Generate Names for Regressor Combinations}
\usage{
get_name_reg_combs(reg_combs)
}
\arguments{
\item{reg_combs}{A nested list of regressor combinations, where each sublist
corresponds to a specific combination size (e.g., all univariate
combinations, all bivariate combinations, etc.), and contains vectors of
regressor names that form each specific combination.}
}
\value{
A list with two elements:
  - `top_lvl_nms`: A vector of names indicating the level of combination
  (univariate, bivariate, etc.) corresponding to the input list's structure.
  - `sub_lvl_nms`: A nested list where each sublist contains names for each
    combination within a level, constructed by concatenating the names of
    regressors in the combination, separated by underscores.
}
\description{
This function constructs meaningful names for each combination of regressors
provided. It assigns names based on the combination's complexity (e.g.,
univariate, bivariate) and concatenates individual regressor names within
each combination for easy identification.
}
\details{
The function is designed to support the identification and tracking of
regressor combinations throughout the data processing and analysis pipeline,
particularly in the context of generating and analyzing multiple regression
models with varying sets of predictors.
}
