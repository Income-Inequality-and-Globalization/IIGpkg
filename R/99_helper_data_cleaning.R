#' Process Covariance Array
#'
#' This function processes a 3D array representing country-period specific 
#' covariance matrices. It symmetrizes the matrices and optionally converts 
#' them into diagonal matrices.
#'
#' @param VCOV_array Initial 3D array of covariance matrices.
#' @param npara Number of parameters (default: 3).
#' @param N Number of countries.
#' @param TT Total number of years.
#' @param Diag_VCOV Boolean flag to convert matrices into diagonal form (default: FALSE).
#' @return A processed 3D array of covariance matrices.
#' @export
process_covariance_array <- function(VCOV_array, npara, N, TT, Diag_VCOV = FALSE) {
  # Symmetrize the covariance matrices
  VCOV_array_pd <- array(apply(VCOV_array, 3,
                               function(x) {
                                 x/2 + t(x)/2
                                }
                               ),
                         c(npara, npara, N * TT))
  # Convert to diagonal matrices if Diag_VCOV is TRUE
  if (Diag_VCOV) {
    VCOV_array_pd <- array(apply(VCOV_array_pd, 3, 
                                 function(x) {
                                   diag(diag(x))
                                  }
                                ),
                           c(npara, npara, N * TT))
  }
  return(VCOV_array_pd)
}

#' Group Data Attributes
#'
#' This function calculates unique years and countries from a dataset and 
#' creates a matrix of country-parameter combinations.
#'
#' @param data_gmm A data frame containing at least the columns 'year' and 'country'.
#' @return A list containing years, total number of years (TT), countries, 
#'         total number of countries (N), and a name matrix (nameMat).
#' @export
get_data_meta_attributes <- function(data_gmm) {
  # Unique years and their count
  years <- unique(data_gmm$year)
  TT <- length(years)

  # Unique countries and their count
  countries <- unique(data_gmm$country)
  N <- length(countries)

  # Creating a matrix for country-parameter combinations
  npara <- 3
  nameMat <- cbind(rep(countries, each = npara), rep(c("a", "q", "mu"), N))

  # Returning the results as a list
  return(list(years = years,
              TT = TT,
              npara = npara,
              countries = countries,
              N = N,
              nameMat = nameMat))
}

#' Remove Specific Year Data and Save
#'
#' This function reads a dataset, filters out a specified year (default 2021),
#' and saves the filtered dataset as an RDS file.
#'
#' @param pth_base_data Base path for data files (default: "./data/input/data-sm").
#' @param data_file Name of the data file (default: "data_coef_covariates.txt").
#' @param save_file Name of the file to save the filtered data (default: "data_coef_covariates_2020.rds").
#' @param year_to_remove Year to be removed from the dataset (default: 2021).
#' @return None
#' @export
remove_year_and_save <- function(pth_base_data = "./data/input/data-sm",
                                 data_file = "data_coef_covariates.txt",
                                 save_file = "data_coef_covariates_2020.rds",
                                 year_to_remove = 2021) {
  pth_data_covr <- file.path(pth_base_data, data_file)
  results_GMM <- tibble::tibble(read.table(pth_data_covr, header = TRUE))
  GMM_by_year <- results_GMM %>% 
                 dplyr::arrange(year) %>% 
                 dplyr::filter(year != year_to_remove)
  saveRDS(GMM_by_year, file = file.path(pth_base_data, save_file))
}
#' Center values by the first non-NA observation
#'
#' This function centers the values in a vector by subtracting the first non-NA 
#' observation from all elements.
#' @param x A numeric vector.
#' @return A numeric vector with values centered around the first non-NA observation.
#'
#' @export
firstObs_center <- function(x){
  nonNA_index <- which(!is.na(x))[1]
  x - x[nonNA_index]
}

#' Retrieve the first non-NA observation
#'
#' This function returns the first non-NA observation from a numeric vector.
#' @param x A numeric vector.
#' @return The first non-NA value in the vector.
#'
#' @export
firstObs_center_values <- function(x){
  x[which(!is.na(x))[1]]
}

#' Standardize a numeric vector
#'
#' This function standardizes a numeric vector by dividing each element by the 
#' standard deviation of the vector, ignoring NA values.
#' @param x A numeric vector.
#' @return A numeric vector with each element standardized.
#'
#' @export
standardize <- function(x){
  x/sd(x, na.rm = TRUE)
}

#' Compute the standard deviation of a numeric vector
#'
#' This function calculates the standard deviation of a numeric vector, 
#' excluding NA values.
#' @param x A numeric vector.
#' @return The standard deviation of the vector, excluding NA values.
#'
#' @export
standardized_values <- function(x){
  sd(x, na.rm = TRUE)
}
