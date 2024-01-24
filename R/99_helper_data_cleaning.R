#' Center values by the first non-NA observation
#'
#' This function centers the values in a vector by subtracting the first non-NA 
#' observation from all elements.
#' @param x A numeric vector.
#' @return A numeric vector with values centered around the first non-NA observation.
#' @examples
#' firstObs_center(c(NA, 2, 3, 4))
firstObs_center <- function(x){
  nonNA_index <- which(!is.na(x))[1]
  x - x[nonNA_index]
}

#' Retrieve the first non-NA observation
#'
#' This function returns the first non-NA observation from a numeric vector.
#' @param x A numeric vector.
#' @return The first non-NA value in the vector.
#' @examples
#' firstObs_center_values(c(NA, 2, 3, 4))
firstObs_center_values <- function(x){
  x[which(!is.na(x))[1]]
}

#' Standardize a numeric vector
#'
#' This function standardizes a numeric vector by dividing each element by the 
#' standard deviation of the vector, ignoring NA values.
#' @param x A numeric vector.
#' @return A numeric vector with each element standardized.
#' @examples
#' standardize(c(1, 2, 3, NA))
standardize <- function(x){
  x/sd(x, na.rm = TRUE)
}

#' Compute the standard deviation of a numeric vector
#'
#' This function calculates the standard deviation of a numeric vector, 
#' excluding NA values.
#' @param x A numeric vector.
#' @return The standard deviation of the vector, excluding NA values.
#' @examples
#' standardized_values(c(1, 2, 3, NA))
standardized_values <- function(x){
  sd(x, na.rm = TRUE)
}
