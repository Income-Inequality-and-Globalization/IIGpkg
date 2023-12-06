#' Adjust for NAs and extend time series for Vhat matrix
#'
#' @param vcov_array_country_pd the true data to adjust and extend
#' @param TT time sereis to extend to
#' @param NN number of countries
#'
#' @return the adjusted dataset
#' @export
adjust_vcov_array_country_pd_to_TT <- function(vcov_array_country_pd, TT, NN) {
  tt_lng <- 21
  id_vcov_na <- which(is.na(vcov_array_country_pd[1,1, ]))
  stopifnot(`Not all entries NA` = all(is.na(vcov_array_country_pd[, , id_vcov_na])))
  
  vcov_array_country_pd_no_na <- vcov_array_country_pd
  for (i in id_vcov_na) {
    vcov_array_country_pd_no_na[, , i] <- vcov_array_country_pd_no_na[, , i - 1]
  }
  stopifnot(`Elimination of NA entries failed` = !any(is.na(vcov_array_country_pd_no_na[, , id_vcov_na])))

  vcov_array_country_pd_extend <- NULL

  tmp_rpl <- ceiling(TT / tt_lng)
  vcov_array_country_pd_extend <- array(replicate(tmp_rpl, vcov_array_country_pd_no_na[, , 1:tt_lng, drop = FALSE]), dim = c(3, 3, tmp_rpl*tt_lng))
  vcov_array_country_pd_extend <- vcov_array_country_pd_extend[, , 1:TT]
  for (nn in 2:NN) {
    vcov_array_country_pd_tkn <- vcov_array_country_pd_no_na[, , ((nn - 1) * tt_lng + 1):(nn * tt_lng)]
    vcov_array_country_pd_tkn <- array(replicate(tmp_rpl, vcov_array_country_pd_tkn),  dim = c(3, 3, tmp_rpl*tt_lng))
    vcov_array_country_pd_tkn <- vcov_array_country_pd_tkn[, , 1:TT]
    vcov_array_country_pd_extend <- abind::abind(
      vcov_array_country_pd_extend,
      vcov_array_country_pd_tkn,
      along = 3
    )
  }
  stopifnot(`Failed to remove and extend NA` = !any(is.na(vcov_array_country_pd_extend)))
  return(vcov_array_country_pd_extend)
}
#' Adjust for zeros and extend time series for regressors
#'
#' @param wregs the true data to adjust and extend
#' @param TT time sereis to extend to
#'
#' @return the adjusted dataset
#' @export
adjust_w_regs_to_TT <- function(wregs, TT) {
  tt_lng <- 21
  wregs_adj_small <- wregs
  wregs_adj_small[, 1] <- wregs_adj_small[, 2]
  wregs_adj_small[28:30, 16:21] <- wregs_adj_small[28:30, 1:6]

  stopifnot(`Failed to eliminate all zeros.` = !any(wregs_adj_small == 0))

  wregs_adj_extend <- wregs_adj_small
  tmp_rpl <- ceiling(TT / tt_lng)
  wregs_adj_extend <- matrix(replicate(tmp_rpl, wregs_adj_extend), nrow = 30)

  wregs_adj_extend <- wregs_adj_extend[, 1:TT]
  rownames(wregs_adj_extend) <- rownames(wregs)
  colnames(wregs_adj_extend) <- paste0("20", formatC(0:(TT - 1), digits = 1, 
                                                     format = "d", flag = "0"))
  return(wregs_adj_extend)
}
