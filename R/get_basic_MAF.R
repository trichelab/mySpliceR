#' Generate a basic MAF-like data.frame from a GRanges object
#'
#' @param gr A GRanges object with metadata columns: ref, alt, and Tumor_Sample_Barcode or sample
#' @param file Optional path to write the resulting data.frame as a tab-delimited MAF file
#'
#' @return A data.frame with the essential MAF columns:
#'   Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, Tumor_Sample_Barcode
#' @export
get_basic_MAF <- function(gr, file = NULL) {
  # Ensure required package is available
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("The 'GenomicRanges' package is required.")
  }

  # Validate GRanges input
  if (!inherits(gr, "GRanges")) {
    stop("Input must be a GRanges object.")
  }

  # Check that ref and alt alleles are present
  required_alleles <- c("ref", "alt")
  if (!all(required_alleles %in% names(mcols(gr)))) {
    stop("GRanges must contain metadata columns: ", paste(required_alleles, collapse = ", "))
  }

  # Determine sample ID column
  sample_col <- if ("Tumor_Sample_Barcode" %in% names(mcols(gr))) {
    "Tumor_Sample_Barcode"
  } else if ("sample" %in% names(mcols(gr))) {
    "sample"
  } else {
    stop("GRanges must contain a 'Tumor_Sample_Barcode' or 'sample' metadata column.")
  }

  # Build the basic MAF data.frame
  maf <- data.frame(
    Chromosome = as.character(seqnames(gr)),
    Start_Position = start(gr),
    Reference_Allele = as.character(mcols(gr)$ref),
    Tumor_Seq_Allele2 = as.character(mcols(gr)$alt),
    Tumor_Sample_Barcode = as.character(mcols(gr)[[sample_col]]),
    stringsAsFactors = FALSE
  )

  # Optionally write to file
  if (!is.null(file)) {
    write.table(maf, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(maf)
}

