% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/99_helper_data_cleaning.R
\name{generate_regressor_combinations}
\alias{generate_regressor_combinations}
\title{Generate Regressor Combinations}
\usage{
generate_regressor_combinations(data_GMM, regs, TT, countries, years)
}
\arguments{
\item{data_GMM}{A data frame containing the regressor data.}

\item{regs}{Vector of regressor names.}

\item{TT}{Total number of years.}

\item{countries}{Vector of unique countries.}

\item{years}{Vector of unique years.}
}
\value{
A list containing two elements: wRegCombsList and nameRegList.
}
\description{
This function generates all possible combinations of given regressors and
creates standardized matrices for each combination.
}
