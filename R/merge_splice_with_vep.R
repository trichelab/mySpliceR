#' Merge mySplice output GRanges with VEP GRanges metadata
#'
#' @param gr_splice A GRanges object returned by mySplice(return_granges = TRUE)
#' @param gr_vep A VEP-annotated GRanges object with names like "chr1:26693521_C/T"
#' @return A GRanges object with an added list-column `vep_info` containing matching VEP metadata rows as GRanges
#' @export
merge_splice_with_vep <- function(gr_splice, gr_vep) {
  if (!inherits(gr_splice, "GRanges") || !inherits(gr_vep, "GRanges")) {
    stop("Both inputs must be GRanges objects")
  }

  # Construct match keys
  key_splice <- paste0(
    as.character(seqnames(gr_splice)), ":",
    start(gr_splice), "_",
    mcols(gr_splice)$ref, "/", mcols(gr_splice)$alt
  )

  key_vep <- names(gr_vep)

  # Get all matches for each splice key and return as GRanges objects
  vep_info <- GenomicRanges::GRangesList(lapply(key_splice, function(k) {
    matched <- which(key_vep == k)
    if (length(matched) > 0) {
      gr_vep[matched]
    } else {
      GenomicRanges::GRanges()
    }
  }))

  # Attach VEP info as a GRangesList column
  mcols(gr_splice)$vep_info <- vep_info

  return(gr_splice)
}

