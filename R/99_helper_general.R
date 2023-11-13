#' Compute Mean and Variance of Inverse Wishart Distribution
#'
#' @param psi Diagonal elements of the Wishart distribution
#' @param nu Degrees of freedom
#' @param p Dimension of the Inverse Wishart (default = 3)
#' @return A list containing the variance and mean of the Inverse Wishart
compute_iw_moments <- function(psi, nu, p = 3) {
  stopifnot(nu > p + 3)  # Ensure nu is greater than p + 3 for valid computation
  mu <- psi/(nu - p - 1)  # Calculate mean
  var_part1 <- (nu - p - 1)^2
  var_part2 <- (nu - p - 3)
  var <- 2 * psi^2/(var_part1 * var_part2)  # Calculate variance
  return(list(var = var, mu = mu))
}
#' Compute Mean and Variance of Inverse Wishart for Multiple Parameters
#'
#' @param psi_taken Vector of diagonal elements of the Wishart distribution
#' @param nu_taken Vector of degrees of freedom
#' @param p Dimension of the Inverse Wishart (default = 3)
#' @return A list containing matrices of means and variances for each psi and nu
#' 
#' @export
#' @examples
#' psi_taken <- c(7, 17, 27, 57)
#' nu_taken <- 7:20
#' iw_moments(psi_taken, nu_taken)
iw_moments <- function(psi_taken, nu_taken, p = 3) {
  out_var <- matrix(nrow = length(psi_taken), ncol = length(nu_taken))
  out_mean <- matrix(nrow = length(psi_taken), ncol = length(nu_taken))
  for (j in seq_along(psi_taken)) {
    out_var[j, ]  <- compute_iw_moments(psi_taken[j], nu_taken, p)[[1]]
    out_mean[j, ] <- compute_iw_moments(psi_taken[j], nu_taken, p)[[2]]
  }
  colnames(out_var) <- paste0("nu_", nu_taken)
  colnames(out_mean) <- paste0("nu_", nu_taken)
  
  rownames(out_var) <- paste0("Psi_", psi_taken)
  rownames(out_mean) <- paste0("Psi_", psi_taken)
  return(list(mean = out_mean, var = out_var))
}