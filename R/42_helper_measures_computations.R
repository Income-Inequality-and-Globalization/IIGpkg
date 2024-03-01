compute_ki_upper_lower <- function(ki_prob) {
  ki_upper <- 1 - (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  return(c(ki_upper, ki_lower))
}
#' Thin Monte Carlo Samples
#'
#' This function applies thinning to Monte Carlo samples specifically from a
#' GibbsOutputIIG class output. Thinning is performed by selecting every nth
#' sample from the Monte Carlo output, where n is defined by the thinning
#' frequency. This method can reduce autocorrelation in the samples and decrease
#' the size of the dataset for analysis.
#'
#' @param x An object of class GibbsOutputIIG, typically the output from a
#'   Monte Carlo sampler. It should contain multidimensional arrays representing
#'   sampled parameter vectors and other relevant matrices.
#' @param thin_freq An integer value representing the thinning frequency, i.e.,
#'   the interval at which samples should be retained. For example, a
#'   thin_freq of 10 means every 10th sample is retained.
#'
#' @return Returns a new GibbsOutputIIG object with thinned samples. The
#'   structure and variables of the input object are preserved, but each
#'   variable is thinned according to the specified frequency.
#'
#' @export
compute_thin_mcmc <- function(x, thin_freq = NULL) {
  if (is.null(thin_freq)) stop("Arg. 'thin_freq' missing or NULL; must be numeric")
  thin_freq_tkn <- floor(thin_freq)
  stopifnot(`Freq. must be an integer` = thin_freq_tkn  == thin_freq)
  validate_GibbsOutputIIG(x)

  MM_MAX_seq <- dim(x$f)[3]
  thin_freq_seq <- seq(from = 1, to = MM_MAX_seq, by  = thin_freq)
  x$f <- x$f[, , thin_freq_seq]
  x$B <- x$B[, , thin_freq_seq]
  x$D <- x$D[, , thin_freq_seq]
  x$A <- x$A[, , thin_freq_seq]
  x$V <- x$V[, thin_freq_seq]
  x$BD0STORE    <- x$BD0STORE[, thin_freq_seq]
  x$Omega0STORE <- x$Omega0STORE[, thin_freq_seq]

  return(new_GibbsOutputIIG(x))
}
# check_unique_vals <- function(array_to_check, dimension) {
#   dim_tmp <- dim(array_to_check)
#
# }
quantile <- function(x, probs) {
  ki_int <- compute_ki_upper_lower(ki_prob = 0.8)
  ki_upp <- ki_int[1]
  ki_low <- ki_int[2]
  stats::quantile(x, probs = c(ki_upp, ki_low))
}