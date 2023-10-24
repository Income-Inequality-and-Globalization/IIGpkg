#' Cleans directory leaving it only with .qmd-files
#'
#' Cleaning means removing everything but the qmd-files.
#'
#' @param dir_to_check character giving the relative path to the directory where
#'    the qmd-files are located
#'
#' @return `NULL`; pure side effect function
#' @export
clean_include_dir <- function(dir_to_check = "./include") {
  # Create temporary directory if it doesn't exist
  tmp_dir <- "./tmp-dir"
  if (!dir.exists(tmp_dir)) {
    dir.create(tmp_dir)
  }

  # List all files and directories in the "include" directory
  files_list <- list.files(dir_to_check, full.names = TRUE)

  # Identify non-qmd files and directories
  files_id <- grep("\\.qmd$", files_list, invert = TRUE, ignore.case = TRUE,)

  for (path in files_list[files_id]) {
    # Generate new path in the temporary directory
    new_path <- file.path(tmp_dir, basename(path))

    # Check if it's a directory
    if (file.info(path)$isdir) {
      # Create a new directory in tmp-dir and copy contents
      if (!dir.exists(new_path)) dir.create(new_path)
      file.copy(path, new_path, recursive = TRUE)
      # Delete the original directory
      unlink(path, recursive = TRUE)
    } else {
      # Move non-qmd files to the temporary directory
      file.rename(path, new_path)
    }

    # Print the new path of moved items
    cat("Moved:", new_path, "\n")
  }

  # Return invisible NULL as the function is for side-effects
  invisible(NULL)
}

#' Cleans the tmp-dir directory
#'
#' Cleaning means removing everything and re-generating an empty dir
#'
#' @param dir_to_check character giving the relative path to the directory where
#'    the qmd-files are located
#'
#' @return `NULL`; pure side effect function
#' @export
clean_tmp_dir <- function(dir_to_check = "./tmp-dir") {
  if (dir.exists(dir_to_check)) unlink(dir_to_check, recursive = TRUE)
  dir.create(dir_to_check)
  return(invisible(NULL))
}
