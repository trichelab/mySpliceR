#' Convert a Gencode GTF file to a BED-format exon junction DataFrame or file
#'
#' @param gtf_path Path to GTF file (required)
#' @param output_bed_path Optional path to save the BED file. If not provided, returns a data.frame
#' @return A data.frame of BED-format exon junctions if `output_bed_path` is not given.
#' @export
#' @examples
#' df <- parse_gtf_to_bed("gencode.v36.annotation.gtf")
#' parse_gtf_to_bed("gencode.v36.annotation.gtf", "junctions.bed")

parse_gtf_to_bed <- function(gtf_path, output_bed_path = NULL) {
  py_file <- system.file("python/parse_gtf_to_junction_bed.py", package = "mySpliceR")
  if (py_file == "") stop("Python script not found in package.")

  reticulate::source_python(py_file)

  if (is.null(output_bed_path)) {
    return(parse_gtf_to_junction_bed(gtf_path))
  } else {
    parse_gtf_to_junction_bed(gtf_path, output_bed_path)
    return(invisible(NULL))
  }
}

