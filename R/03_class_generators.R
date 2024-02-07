#' Class generator for Gibbs output
#'
#' @param x output object (a named list) as returned via
#'   [IIGpkg::Gibbs2_SM_SA_sampler] and [IIGpkg::GibbsSSM_2]
#'
#' @return an object of type `GibbsOutputIIG` which is a named list; see the 
#'   return value of [IIGpkg::GibbsSSM_2] for details.
#' @export
new_GibbsOutputIIG <- function(x = vector("list")) {
  validate_GibbsOutputIIG(x)
  structure(x, class = "GibbsOutputIIG")
}
validate_GibbsOutputIIG <- function(x) {
  stopifnot(is.list(x))
  stopifnot(`List names malforamtted` = (all(names(x) %in% c("f", "B", "D",
                                        "A", "V", "u",
                                        "BD0STORE", "Omega0STORE",
                                        "blockCount",
                                        "errorMsg",
                                        "initials"))))
  return(invisible(NULL))
}
#' Subset Gibbs output with smaller than total Monte Carlo size
#'
#' Some of the Gibbs outputs are very large and it is not always necessary to 
#' use all of the Monte Carlo draws. Hence, by passing a numeric vector sequence
#' number `MM_seq` in the range of total Monte Carlo draws, we can subset the
#' output for all parameters thus making the object smaller in size and related
#' computations (Bayesian averages or confidence bands) faster to compute.
#'
#' @inheritParams new_GibbsOutputIIG
#' @param MM_seq numeric vector; vector of Monte Carlo indices to keep which 
#'   should be in the interval between 1 and the total number of draws (of the
#'   container objects stored in  `x`)
#' 
#' @return same as `x` but with the Monte Carlo draws taken from the range in
#'    `MM_seq`
#' @export
subset_GibbsOutputIIG <- function(x, MM_seq = NULL) {
  validate_GibbsOutputIIG(x)
  if (is.null(MM_seq)) stop("Arg. 'MM' missing or NULL; must be numeric")
  MM_MAX_seq <- seq_len(dim(x$f)[3])
  stopifnot(`Arg. 'MM_seq' sequence not in range` = all(MM_seq %in% MM_MAX_seq))

  x   <- x
  x$f <- x$f[, , MM_seq]
  x$B <- x$B[, , MM_seq]
  x$D <- x$D[, , MM_seq]
  x$A <- x$A[, , MM_seq]
  x$V <- x$V[, MM_seq]
  x$BD0STORE    <- x$BD0STORE[, MM_seq]
  x$Omega0STORE <- x$Omega0STORE[, MM_seq]

  return(new_GibbsOutputIIG(x))
}
#' Summary of Gibbs output
#' 
#' A table with true values (if they exist), Monte Carlo means, confidence
#' bands, and an indicator if the true value lies in the CI.
#'
#' @inheritParams new_GibbsOutputIIG
#' @param true_vals a list of object containing the true values
#' @param options a list of options:
#'   * `nrows_out`: number of rows to print for summary
#'   * `roundoff`: integer giving the digits to round values in summary
#'   * `ci_val`: value of confidence band; default is 1.96 for 95\%
#'
#' @return pure side effect functions but returns invisibly a list of the 
#'   different outputs that are printed to the screen
#' @export
summary.GibbsOutputIIG <- function(x, true_vals = NULL, burnin = NULL,
                                   options = list(nrows_out = 9999,
                                                  sig_digits = 4,
                                                  roundoff = 2,
                                                  ci_val = 1.96)) {
  old <- options(pillar.sigfig = options$sig_digits)
  on.exit(options(old))
  validate_GibbsOutputIIG(x)

  true_B_means <- true_vals$B_means

  num_jnt_fac <- x$initials$njointfac
  num_all_facs <- ncol(x$B)
  # num_regs    <- x$initials$nreg
  EXIST_JNT_FAC <- num_jnt_fac != 0
  
  num_mcmc <- x$initials$itermax
  if (is.null(burnin)) {
    num_burn <- floor(num_mcmc / 2)
  } else {
    num_burn <- burnin
  }

  out <- x$B[, (num_jnt_fac + 1):num_all_facs, (num_burn + 1):num_mcmc]
  mean_il <- diag(apply(out, c(1,2), mean))
  if (!is.null(true_vals)) {
      CIs <- apply(out, 3, diag)
      CIs <- apply(CIs, 1, quantile, probs = c(0.025, 0.975))
      CI_lower <- CIs[1, ]
      CI_upper <- CIs[2, ]
      true_means <- diag(true_B_means[, (num_jnt_fac + 1):num_all_facs])
      contained <- (true_means <= CI_upper) & (true_means >= CI_lower)
  } else {
      true_means <- rep(NA_real_, times = length(mean_il))
      contained <- true_means
  }
  out_il <- tibble::tibble(true = round(true_means, digits = options$roundoff),
                           `post. mean` = round(mean_il, digits = options$roundoff),
                           CI_lower = round(CI_lower, digits = options$roundoff),
                           CI_upper = round(CI_upper, digits = options$roundoff),
                           `if in CI -> 1` = round(contained, digits = options$roundoff))
  message(crayon::yellow("Summary for idiosyncratic loadings: "))
  print(out_il, n = options$nrows_out)
    
  if (EXIST_JNT_FAC) {
    out2 <- x$B[, 1:num_jnt_fac, (num_burn + 1):num_mcmc]
    mean_cl <- apply(out2, 1, mean)
    if (!is.null(true_vals)) {
      CIs <- apply(out2, 1, quantile, probs = c(0.025, 0.975))
      CI_lower <- CIs[1, ]
      CI_upper <- CIs[2, ]
      true_means <- true_B_means[, 1:num_jnt_fac]
      contained <- (true_means <= CI_upper) & (true_means >= CI_lower)
    } else {
      true_means <- rep(NA_real_, times = length(mean_il))
      contained <- true_means
    }
    out_cl <- tibble::tibble(true = round(true_means, digits = options$roundoff),
                             `post. mean` = round(mean_cl, digits = options$roundoff),
                             CI_lower = round(CI_lower, digits = options$roundoff),
                             CI_upper = round(CI_upper, digits = options$roundoff),
                             `if in CI -> 1` = round(contained, digits = options$roundoff))
    message(crayon::yellow("Summary for joint loadings: "))
    print(out_cl, n = options$nrows_out )
  }
  options(old)
  return(list(summary_loadings_idiosyncratic = out_il,
              summary_loadings_joint = out_cl))
}
