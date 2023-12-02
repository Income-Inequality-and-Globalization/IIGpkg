#' Generate a list with relevant true (simulation) values
#'
#' @param pth_to_true_vals pth where the `.rds` file is stored
#' @param ... the relevant values to set, see function body
#'
#' @return invisibly saving an `.rds` file
#' @export
#'
#' @examples 
#' pth_to_true_vals("./kf_test2.RData")
#' set_true_vals(pth_to_true_vals,
#'               VhatArray_A,
#'               A_means,
#'               B_means,
#'               D_means,
#'               yObs,
#'               w_regs,
#'               f_simul)
set_true_vals <- function(pth_to_true_vals,
                          ...) {
  names_to_set <- list(...)
  tmp_names <- c("VhatArray_A",
                 "A_means",
                 "B_means",
                 "D_means",
                 "y_used",
                 "w_regs_used",
                 "f_simul")
  out <- vector("list", length = length(tmp_names))
  
  out$VhatArray_A <- names_to_set[[1]]
  out$A_means     <- names_to_set[[2]]
  out$B_means     <- names_to_set[[3]]
  out$D_means     <- names_to_set[[4]]
  out$y_used      <- names_to_set[[5]]
  out$w_regs_used <- names_to_set[[6]]
  out$f_simul     <- names_to_set[[7]]
  
  all_true_vals <- out
  save(all_true_vals, file = pth_to_true_vals)
  return(invisible(all_true_vals))
}