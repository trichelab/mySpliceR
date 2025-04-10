#' Extract all read names from a mySplice GRanges object
#'
#' @param gr GRanges object with a 'supporting_reads' column
#' @return A named list of character vectors, one per variant
#' @export
get_all_read_names <- function(gr) {
  if (!inherits(gr, "GRanges")) stop("gr must be a GRanges object.")
  if (!"supporting_reads" %in% colnames(mcols(gr))) stop("'supporting_reads' column not found.")

  out <- lapply(seq_along(gr), function(i) {
    reads_df <- gr$supporting_reads[[i]]
    if (is.null(reads_df) || nrow(reads_df) == 0) return(character(0))
    if (!"read_name" %in% colnames(reads_df)) stop("'read_name' column not found in supporting_reads.")
    unique(reads_df$read_name)
  })

  names(out) <- if (!is.null(names(gr))) names(gr) else as.character(seq_along(gr))
  return(out)
}

