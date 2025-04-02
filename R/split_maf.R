#' Split a MAF file into smaller chunks
#'
#' This function splits a MAF (Mutation Annotation Format) file into approximately equal-sized chunks,
#' while preserving the header line in each output file. The splitting is handled by a Python script
#' bundled in the `inst/python` directory of the package.
#'
#' @param file_in Path to the input MAF file.
#' @param dir_out Directory where the split files will be written. It will be created if it doesn't exist.
#' @param max_num Maximum number of output files to create.
#'
#' @return Invisibly returns the system exit code (0 if successful).
#' @examples
#' \dontrun{
#'   split_maf("inst/extdata/example.maf", "inst/extdata/splits", 5)
#' }
#' @export
split_maf <- function(file_in, dir_out, max_num) {
  script_path <- system.file("python", "split_maf.py", package = "mySpliceR")

  if (!file.exists(script_path)) {
    stop("Could not find split_maf.py in the package. Please check installation.")
  }

  # Ensure output directory exists
  if (!dir.exists(dir_out)) {
    dir.create(dir_out, recursive = TRUE)
  }

  # Build command
  cmd <- sprintf("python3 %s %s %s %s",
                 shQuote(script_path),
                 shQuote(file_in),
                 shQuote(dir_out),
                 shQuote(as.character(max_num)))

  # Run the command
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("split_maf.py exited with non-zero status: ", exit_code)
  }

  invisible(exit_code)
}

