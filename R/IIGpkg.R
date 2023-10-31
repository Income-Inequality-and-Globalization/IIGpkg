#' IIGpkg: Package for replicating the IIG project results
#'
#' This package provides facilities for replicating analysis of the IIG
#' project. Regarding the compilation of docs, see `Startup`.
#'
#' @section Startup:
#'     To quickly start a new project run:
#' \itemize{
#'  \item{1.}{Copy paste this package to the new location and open `.Rproj` to
#'  get inside the package project with `getwd()` set correctly.}
#'  \item{2.}{Run `library(quartoWebsiteBlueprint)`
#'     and `setup_new_webpage(project_name = new_project_name)`, leaving the
#'     first argument `path` to the default `getwd()` (due to 1.).
#'     This sets a new name for the project that is different from this default
#'     name "quartoWebsiteBlueprint".}
#'  \item{3.}{ Run `quarto render` in the command line. Then make a new commit.}
#'  \item{4.}{ Add github repository, follow the outlined rules to add this
#'  remote to the local version (rename branch-name to `main` and push with
#'  `--set-upstream`).}
#' }
#'
#' @docType package
#' @name IIGpkg
NULL
