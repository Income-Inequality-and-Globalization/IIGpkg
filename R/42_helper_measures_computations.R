compute_ki_upper_lower <- function(ki_prob) {
  ki_upper <- 1 - (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  return(c(ki_upper, ki_lower))
}
get_id_D_uncertain <- function(all_names_Dcol, name_reg) {
  which(grepl(name_reg, all_names_Dcol))
}
get_D_uncertain <- function(D, D_avg, mm, id_cu) {
  D_avg[, id_cu, mm] <- D[, id_cu, mm]
  return(D_avg[, , mm])
}