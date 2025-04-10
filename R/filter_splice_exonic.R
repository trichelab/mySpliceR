#' Filter mySplice GRanges to retain only variants with at least one exon annotation in VEP
#'
#' @param gr_splice A GRanges object returned by merge_splice_with_vep(), with a vep_info GRangesList column
#' @param exon_keywords A character vector of terms considered as exon annotations (default includes \"exon\")
#' @return A filtered GRanges object
#' @export
filter_splice_exonic <- function(gr_splice, exon_keywords = c("exon", "exonic", "protein_coding", "missense", "nonsense")) {
  if (!"vep_info" %in% colnames(mcols(gr_splice))) {
    stop("gr_splice must contain a 'vep_info' column (GRangesList)")
  }

  # Check each variant's vep_info for an exon-related consequence
  has_exon_hit <- vapply(gr_splice$vep_info, function(gr) {
    if (length(gr) == 0) return(FALSE)
    if (!"Consequence" %in% colnames(mcols(gr))) return(FALSE)
    any(grepl(paste(exon_keywords, collapse = "|"), mcols(gr)$Consequence, ignore.case = TRUE))
  }, logical(1))

  # Filter
  gr_splice[has_exon_hit]
}
