# zzz.R — lifecycle hooks for mySpliceR
.onLoad <- function(libname, pkgname) {
  missing <- c()

  for (pkg in c("pandas", "pysam", "shutil")) {
    if (!reticulate::py_module_available(pkg)) {
      missing <- c(missing, pkg)
    }
  }

  if (length(missing) > 0) {
    packageStartupMessage(
      "The following Python packages are missing: ", paste(missing, collapse = ", "), "\n",
      "Use mySpliceR::install_python_dependencies() to install them."
    )
  }
}


