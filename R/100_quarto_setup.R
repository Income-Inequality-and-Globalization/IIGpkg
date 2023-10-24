#' Setup helper for webpage
#'
#' A wrapper around [changer::changer()] that changes the package name and `git`
#' repository of this blueprint. Helps to quickly publish the new webpage
#' without having to rename all files.
#'
#' @param project_name the new project or package name
#' @param path to the blueprint quarto package
#'
#' @return pure side effect function; returns invisibly
#' @export
setup_new_webpage <- function(project_name, path = NULL) {
  if (is.null(path)) path <- getwd()
  pth_to_project_old <- path
  pth_to_project_new <- file.path(dirname(pth_to_project_old), project_name)
  cat(crayon::green("Change from: "), crayon::yellow(pth_to_project_old), "\n")
  cat(crayon::green("to new project: "), crayon::yellow(pth_to_project_new), "\n")
  changer::changer(path = path,
                   new_name = project_name,
                   ask = FALSE)

  name_main_file <- paste0(file.path(getwd(), paste0(project_name, ".qmd")))
  file.rename(paste0(file.path(getwd(), "quartoWebsiteBlueprint.qmd")),
              name_main_file)

  unlink(pth_to_project_old, recursive = TRUE)

  stopifnot(getwd() == pth_to_project_new)
  tmp_git_remote <- gert::git_remote_list()
  unlink(".git", recursive = TRUE)
  gert::git_init()
  quarto::quarto_render(name_main_file)
  gert::git_add(".")
  gert::git_commit("initial commit")
  gert::git_branch_move("master", "main")
  gert::git_remote_add(url = tmp_git_remote$url, tmp_git_remote$name)

  rstudioapi::openProject(pth_to_project_new, newSession = TRUE)
}
#' Clean up artifact directories/files
#'
#' After first creation/renaming via [setup_new_webpage()] there remains the old
#' artifact directory '.../quartoWebsiteBlueprint'. It arises probably due to
#' some `RStudio` internal processes and should be deleted. This function does
#' this safely.
#'
#' @param path a character string; the default value should be fine
#'
#' @return pure side effect function; returns invisibly
#' @export
clean_artifacts <- function(path = getwd()) {
  pth_delete <- file.path(dirname(path), "quartoWebsiteBlueprint")
  if (isFALSE(dir.exists(pth_delete))) {
    stop("Aborting cleanup of artifacts: could not find artifact dir.")
  }
  message(
    paste0(
      "The directory contains the following files:\n",
      list.files(pth_delete)
    )
  )
  msg <- paste0("Delete directory: ", pth_delete, " ? ")
  YES <- utils::askYesNo(msg)
  if (isTRUE(YES)) {
    unlink(pth_delete, recursive = TRUE)
  } else {
    message("Cleanup aborted by user.")
  }
  return(invisible(path))
}
