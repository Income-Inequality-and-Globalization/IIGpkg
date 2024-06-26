Sys.setenv(VSCODE_DEBUG = TRUE)
env <- Sys.getenv()
envnames <- names(env)
rnames <- envnames[startsWith(envnames, "R_")]
cached_names <- rnames
ld_lib_path <- Sys.getenv("LD_LIBRARY_PATH")
if (ld_lib_path != "") {
  cached_names <- c("LD_LIBRARY_PATH", rnames)
}
if (Sys.getenv("VSCODE_DEBUG") == "TRUE") {
  cached_names <- c("VSCODE_DEBUG", cached_names)
}
writeLines(paste0(cached_names, "=", env[cached_names]), ".vscode/.env")