#' Helper to pre-compute adjustment because of missing values in a KF
#'
#' @param yObs Y measurements/observations with possible missing values
#' @param DwReg DxwReg matrices to be adjusted to the missings in Y
#' @param R R matrix to be adjusted to the missings in Y
#' @param C C matrix to be adjusted to the missings in Y
#' @param TT time series length
#' @param dimY length of measurement component
#'
#' @return a list of four elements being lists of length `TT` with adjusted 
#'    matrices (for C, R, yObs, and `D` times `wReg`)
#' @export
get_ffbs_missing_obs_adj <- function(yObs, DwReg, R, C, TT, dimY) {
  out_CAdj     <- list(TT)
  out_RAdj     <- list(TT)
  out_yObsAdj  <- list(TT)
  out_DwRegAdj <- list(TT)
  for (t in 1:TT) {
    AllObsMissing <- FALSE
    CAdj <- C
    if (is.na(dim(R)[3])) {
      RAdj <- R
    } else {
      RAdj <- R[, , t]
    }
    yObsAdj <- yObs[, t]
    DwRegAdj <- DwReg[, t]
    MissingObs <- which(is.na(yObsAdj))
    
    if (length(MissingObs) != 0) {
      if (length(MissingObs) == dimY) {
        AllObsMissing <- TRUE
      } else {
        W <- diag(dimY)
        W <- W[-MissingObs, ]
        
        yObsAdj[MissingObs] <- 0
        CAdj[MissingObs, ] <- 0
        RAdj[MissingObs, ] <- 0
        RAdj[ ,MissingObs] <- 0
        
        out_CAdj[[t]] <- W %*% CAdj
        out_RAdj[[t]] <- W %*% RAdj %*% t(W)
        out_yObsAdj[[t]] <- W %*% yObsAdj
        out_DwRegAdj[[t]] <- W %*% DwRegAdj
      }
    } else {
      out_CAdj[[t]] <- CAdj
      out_RAdj[[t]] <- RAdj
      out_yObsAdj[[t]] <- yObsAdj
      out_DwRegAdj[[t]] <- DwRegAdj
    }
  }
  return(list(CAdj = out_CAdj, RAdj = out_RAdj,
              yObsAdj = out_yObsAdj, DwRegAdj = out_DwRegAdj))
}
adjust_c <- function(CAdj, dimY, MissingObs) {
  CAdj2 <- CAdj
  for (t in 1:dimY) {
    for (j in MissingObs) {
      CAdj2[j, t] <- 0
    }
  }
  return(CAdj2)
}
adjust_r <- function(RAdj, MissingObs) {
  RAdj2 <- RAdj
  for (j in MissingObs) {
    for (k in MissingObs) {
      RAdj2[j, k] <- 0
    }
  }
  return(RAdj2)
}
