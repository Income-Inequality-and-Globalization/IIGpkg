get_subset_par <- function(out_post_all,
                           par_name,
                           SUBS_ROWS = TRUE,
                           SUBS_COLS = FALSE,
                           id_keep_cols = NULL,
                           id_keep_rows = NULL) {
  stopifnot(`Invalid argument 'par_name'.` = par_name %in% c("a", "q", "mu"))
  if (SUBS_ROWS) {
    row_names_par <- rownames(out_post_all)
    id_taken <- find_id_par_name_all(par_name, row_names_par)
    if (!is.null(id_keep_rows)) id_taken <- sort(c(id_taken, id_keep_rows))
    out_post_all <- out_post_all[id_taken, , ]
  }
  if (SUBS_COLS) {
    col_names_par <- colnames(out_post_all)
    id_taken <- find_id_par_name_all(par_name, col_names_par)
    if (!is.null(id_keep_cols)) id_taken <- sort(c(id_taken, id_keep_cols))
    out_post_all <- out_post_all[, id_taken, ]
  }
  return(out_post_all)
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