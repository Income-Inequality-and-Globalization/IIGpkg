set_simulate_factors <- function(fPost) {
  if (missing(fPost)) return(TRUE)
  FALSE
}
set_reg_specification <- function(wRegSpec) {
  if (missing(wRegSpec)) wRegSpec <- 0
  wRegSpec
}