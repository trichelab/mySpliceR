#' Install Python dependencies for mySpliceR
#' @param method Installation method ("auto", "virtualenv", or "conda")
#' @export
install_python_dependencies <- function(method = "auto") {
  pkgs <- c("pandas", "pysam")
  # shutil and os are built-in, no need to install
  reticulate::py_install(packages = pkgs, method = method)
}
