get_par_post_estim <- function(out_post_all, par_name) {
  row_names_par <- rownames(out_post_all)
  id_taken <- find_id_par_name_all(par_name, row_names_par)
  out_post_all[id_taken, , ]
}
find_id_par_name_all <- function(par_name, names_to_search) {
  tmp_rgx <- get_regex_par(par_name)
  which(grepl(tmp_rgx, names_to_search))
}
get_regex_par <- function(par_name) {
  if (par_name == "a") return("a")
  if (par_name == "q") return("q")
  if (par_name == "mu") return("mu")
}