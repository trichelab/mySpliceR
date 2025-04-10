#' Extract read names for a single variant row in mySplice GRanges
#'
#' @param gr GRanges object with a 'supporting_reads' column
#' @param row_index Integer index of the row (variant)
#' @return A character vector of read names
#' @export
get_read_names_for_row <- function(gr, row_index) {
  if (!inherits(gr, "GRanges")) stop("gr must be a GRanges object.")
  if (!"supporting_reads" %in% colnames(mcols(gr))) stop("'supporting_reads' column not found.")
  if (row_index < 1 || row_index > length(gr)) stop("row_index out of bounds.")

  reads_df <- gr$supporting_reads[[row_index]]
  if (is.null(reads_df) || nrow(reads_df) == 0) return(character(0))
  if (!"read_name" %in% colnames(reads_df)) stop("'read_name' column not found in supporting_reads.")

  return(unique(reads_df$read_name))
}

