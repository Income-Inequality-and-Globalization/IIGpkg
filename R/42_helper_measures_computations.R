compute_ki_upper_lower <- function(ki_prob) {
  ki_upper <- 1 - (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  ki_lower <- (1 - ki_prob) / 2
  return(c(ki_upper, ki_lower))
}