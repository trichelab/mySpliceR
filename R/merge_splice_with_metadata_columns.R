#' Merge metadata columns from a GRanges to the GRanges from mySplice (1-to-1 match)
#'
#' @param gr_splice A GRanges object from mySplice(return_granges = TRUE)
#' @param gr_meta A GRanges with columns: sample, ref, alt, and metadata to attach
#' @return A GRanges object with metadata columns from gr_meta added
#' @export
merge_splice_with_metadata_columns <- function(gr_splice, gr_meta) {
  if (!inherits(gr_splice, "GRanges") || !inherits(gr_meta, "GRanges")) {
    stop("Both inputs must be GRanges objects")
  }

  # Construct keys
  key_splice <- paste0(
    gr_splice$sample, "_",
    as.character(seqnames(gr_splice)), "_",
    start(gr_splice), "_",
    gr_splice$ref, "_", gr_splice$alt
  )
  key_meta <- paste0(
    gr_meta$sample, "_",
    as.character(seqnames(gr_meta)), "_",
    start(gr_meta), "_",
    gr_meta$ref, "_", gr_meta$alt
  )

  # Match indices
  matched_idx <- match(key_splice, key_meta)

  # Extract matched metadata
  matched_meta <- mcols(gr_meta)[matched_idx, , drop = FALSE]

  # Optional: check for any NA matches
  if (any(is.na(matched_idx))) {
    warning("Some variants in gr_splice did not match gr_meta. NAs introduced.")
  }

  # Combine metadata columns into gr_splice
  mcols(gr_splice) <- cbind(mcols(gr_splice), matched_meta)

  return(gr_splice)
}

