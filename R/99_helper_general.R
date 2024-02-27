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
#' Compute Mean and Variance of Inverse Gamma Distribution
#'
#' @param alpha Shape parameter of the Inverse Gamma distribution
#' @param beta Scale parameter of the Inverse Gamma distribution
#' @return A list containing the variance and mean of the Inverse Gamma
compute_ig_moments <- function(alpha, beta) {
  stopifnot(alpha > 2)  # Ensure alpha is greater than 2 for valid computation
  mu <- beta / (alpha - 1)  # Calculate mean
  var <- beta^2 / ((alpha - 1)^2 * (alpha - 2))  # Calculate variance
  return(list(mean = mu, variance = var))
}
#' Compute Mean and Variance of Inverse Gamma for Multiple Parameters
#'
#' @param alpha_taken Vector of shape parameters of the Inverse Gamma 
#'    distribution
#' @param beta_taken Vector of scale parameters
#' @return A list containing matrices of means and variances for each alpha and
#'    beta
#'
#' @export
#' @examples
#' alpha_taken <- c(3, 5, 7, 10)
#' beta_taken <- c(2, 4, 6, 8)
#' ig_moments(alpha_taken, beta_taken)
ig_moments <- function(alpha_taken, beta_taken) {
  out_var <- matrix(nrow = length(alpha_taken), ncol = length(beta_taken))
  out_mean <- matrix(nrow = length(alpha_taken), ncol = length(beta_taken))
  for (i in seq_along(alpha_taken)) {
    for (j in seq_along(beta_taken)) {
      moments <- compute_ig_moments(alpha_taken[i], beta_taken[j])
      out_mean[i, j] <- moments$mean
      out_var[i, j] <- moments$variance
    }
  }
  colnames(out_var) <- paste0("Beta_", beta_taken)
  colnames(out_mean) <- paste0("Beta_", beta_taken)
  
  rownames(out_var) <- paste0("Alpha_", alpha_taken)
  rownames(out_mean) <- paste0("Alpha_", alpha_taken)
  return(list(mean = out_mean, variance = out_var))
}
progress_any <- function(iter,
                         iter_max,
                         settings = list(digits = 0,
                                         repeat_every = 10)) {
  percentage_progress <- iter / iter_max * 100
  
  chunck_to_print  <- floor(iter_max / settings$repeat_every)
  repeat_every_seq <- seq(from = chunck_to_print,
                          to = iter_max,
                          by = chunck_to_print)
  check_print_progress <- iter %in% repeat_every_seq
  if (check_print_progress) {
    cat(crayon::green("Progress:"),
        crayon::yellow(
          paste0(round(percentage_progress, digits = settings$digits), "%.\n")))
  }
}
