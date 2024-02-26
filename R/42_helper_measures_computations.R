compute_ki_upper_lower <- function(ki_prob) {
  ki_upper <- 1 - (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  return(c(ki_upper, ki_lower))
}
get_col_reg_name <- function(all_names_Dcol, name_reg) {
  which(grepl(name_reg, all_names_Dcol))
}
get_D_uncertain <- function(D, D_avg, id_cu) {
  D_avg[, id_cu] <- D[, id_cu]
  return(D_avg)
}