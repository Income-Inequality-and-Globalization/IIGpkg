#' Generate Regressor Combinations
#'
#' This function generates all possible combinations of given regressors and
#' creates standardized matrices for each combination.
#'
#' @param data_GMM A data frame containing the regressor data.
#' @param regs Vector of regressor names.
#' @param TT Total number of years.
#' @param countries Vector of unique countries.
#' @param years Vector of unique years.
#'
#' @return A list containing two elements: wRegCombsList and nameRegList.
#' @export
generate_regressor_combinations <- function(data_GMM,
                                            regs,
                                            TT,
                                            countries,
                                            years) {
  lregs <- length(regs)
  reg_combs <- lapply(1:lregs, \(x) combn(regs, x, simplify = FALSE))
  ncombs <- sum(sapply(reg_combs, length))

  wRegCombsList <- vector("list", ncombs)
  nameRegList <- vector("list", ncombs)

  ii <- 1
  for (i in 1:length(reg_combs)) {
    for (j in 1:length(reg_combs[[i]])) {
      selected_regressors <- reg_combs[[i]][[j]]
      wReg <- matrix(t(data_GMM %>% 
                         dplyr::select(tidyselect::all_of(selected_regressors))),
                     ncol = TT)
      wReg_centered <- t(apply(wReg, 1, firstObs_center))
      wReg <- t(apply(wReg_centered, 1, standardize))
      wReg[is.na(wReg)] <- 0
      rownames(wReg) <- rep(countries, each = length(selected_regressors))
      colnames(wReg) <- years
      nameRegList[[ii]] <- substr(selected_regressors, 1, 7)
      wRegCombsList[[ii]] <- wReg
      ii <- ii + 1
    }
  }

  return(list(wRegCombsList = wRegCombsList, nameRegList = nameRegList))
}

#' Adjust Covariance Matrix with Logarithmic Transformation
#'
#' This function adjusts a 3D array of covariance matrices using logarithmic 
#' transformations of the observation data. It symmetrizes the adjusted
#' matrices.
#'
#' @param VCOV_array Initial 3D array of covariance matrices.
#' @param npara Number of parameters.
#' @param N Number of countries.
#' @param TT Total number of years.
#' @param yObs_log Logarithmically transformed observation data.
#' @return A 3D array of adjusted and symmetrized covariance matrices.
#' @export
log_adj_VCOV <- function(VCOV_array, npara, N, TT, yObs_log) {
  yObs_inv <- 1 / yObs_log
  ArrayTimeCountry <- array(yObs_inv, c(npara, 1, N * TT))
  MatCountryTime <- ArrayTimeCountry[, , c(sapply(0:(N - 1), \(x) (x + seq(1, (TT - 1) * N + 1, N))))]
  VCOV_log_list <- lapply(
    1:(N * TT),
    \(x) diag(MatCountryTime[, x]) %*% VCOV_array[, , x] %*% diag(MatCountryTime[, x]))
  VCOV_log_array <- array(unlist(VCOV_log_list), dim = c(npara, npara, TT * N))
  array(apply(VCOV_log_array, 3, function(x) {x / 2 + t(x) / 2}), c(npara, npara, N * TT))
}

#' Standardize Covariance Matrices
#'
#' This function standardizes a 3D array of covariance matrices using the 
#' standard deviations of the logarithmically transformed observation data.
#'
#' @param VCOV_array A 3D array of covariance matrices to be standardized.
#' @param TT Total number of years.
#' @param N Number of countries.
#' @param sd_vec_y Vector of standard deviations for each parameter and country.
#' @param npara Number of parameters.
#' @return A 3D array of standardized covariance matrices.
#' @export
standardize_VCOV <- function(VCOV_array, TT, N, sd_vec_y, npara) {
  V_adj_array <- array(0, dim = c(npara, npara, TT * N))
  for (i in 1:N) {
    V_adj <- diag(1 / sd_vec_y[((i - 1) * npara + 1):(i * npara)])
    V_select <- VCOV_array[, , ((i - 1) * TT + 1):(i * TT)]
    adj_array <- array(apply(V_select, 3, \(x) V_adj %*% x %*% V_adj),
                       dim = c(npara, npara, TT))
    V_adj_array[, , ((i - 1) * TT + 1):(i * TT)] <- adj_array
  }
  return(V_adj_array)
}

#' Process Observation Data
#'
#' This function processes observation data from a dataset, performing log
#' transformation, centering, standardization, and calculating standard
#' deviations.
#'
#' @param data_GMM A data frame with columns 'a', 'q', 'est_mean', 'country',
#'   and 'year'.
#' @param npara Number of parameters (default: 3).
#' @param TT Total number of years.
#' @param countries Vector of unique countries.
#' @param years Vector of unique years.
#' @return A list containing various transformations of the observation data.
#'
#' @export
process_observation_data <- function(data_GMM, npara, TT, countries, years) {
  # Create the unadjusted observation matrix
  yObs_unadj <- matrix(t(data_GMM %>% dplyr::select(a, q, est_mean)), ncol = TT)
  rownames(yObs_unadj) <- rep(countries, each = npara)
  colnames(yObs_unadj) <- years

  # Log transformation
  yObs_log <- log(yObs_unadj)

  # Centering and standardization
  yObs_log_centered <- t(apply(yObs_log, 1,
                               firstObs_center))
  yObs_log_centered_values <- t(apply(yObs_log, 1,
                                      firstObs_center_values))
  yObs_log_centered_standardized <- t(apply(yObs_log_centered, 1,
                                            standardize))
  yObs_standardized_values <- t(apply(yObs_log_centered, 1, 
                                      standardized_values))

  # Standard deviation
  sd_yObs_log <- apply(yObs_log, 1, sd, na.rm = TRUE)
  # Return all variables as a list
  return(
    list(
      data = list(
        yObs_unadj = yObs_unadj,
        yObs_log = yObs_log, 
        yObs_log_centered = yObs_log_centered, 
        yObs_log_centered_standardized = yObs_log_centered_standardized),
      data_values = list(
        yObs_log_centered_values = yObs_log_centered_values,
        yObs_standardized_values = yObs_standardized_values, 
        sd_yObs_log = sd_yObs_log)
      )
  )
}

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
#' @param Diag_VCOV Boolean flag to convert matrices into diagonal form
#'   (default: FALSE).
#' @return A processed 3D array of covariance matrices.
#' @export
process_covariance_array <- function(VCOV_array,
                                     npara, N, TT,
                                     Diag_VCOV = FALSE) {
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
#' @param data_gmm A data frame containing at least the columns 'year' and
#'   'country'.
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
#' and saves the filtered dataset as an `.rds` file. It also adds the ginis for
#' countries and years and saves the data to a `.txt`-file
#'
#' @param pth_base_data Base path for data files (default:
#'   "./data/input/data-sm").
#' @param data_file Name of the data file (default: "data_coef_covariates.txt").
#' @param save_file Name of the file to save the filtered data (default:
#'   "data_coef_covariates_2020.rds").
#' @param year_to_remove Year to be removed from the dataset (default: 2021).
#' @return None
#' @export
data_raw_remove_year_and_save <- function(
    pth_base_data = "./data/input/data-sm",
    data_file = "data_coef_covariates.txt",
    save_file_01 = "data_coef_covariates_2020.rds",
    save_file_02 = "data_coef_covariates_with_ginis_2020.csv",
    year_to_remove = 2021) {
  pth_data_covr <- file.path(pth_base_data, data_file)
  results_GMM  <- tibble::tibble(read.table(pth_data_covr, header = TRUE))
  results_gini <- tibble::tibble(read.csv(file.path(pth_base_data,
                                                    "../raw-sources-others",
                                                    "dataset_add_regs.csv")))
  data_full <- dplyr::full_join(results_GMM, results_gini)
  data_GMM  <- results_GMM %>% 
                  dplyr::arrange(year) %>% 
                  dplyr::filter(year != year_to_remove)
  data_full  <- data_full %>% 
                  dplyr::arrange(country, year) %>% 
                  dplyr::filter(year != year_to_remove) %>%
                  dplyr::mutate(gini_from_gmm = compute_gini_sm(a, q) * 100)
  saveRDS(data_GMM, file = file.path(pth_base_data, save_file_01))
  write.csv(data_full, 
            file = file.path(pth_base_data, save_file_02),
            row.names = FALSE)
}
#' Center values by the first non-NA observation
#'
#' This function centers the values in a vector by subtracting the first non-NA 
#' observation from all elements.
#' @param x A numeric vector.
#' @return A numeric vector with values centered around the first non-NA
#'   observation.
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
