#' Detect mis-splicing events based on somatic mutations
#'
#' @param maf_file Path to the input MAF file
#' @param bam_list Path to the TSV file listing sample ID and BAM paths
#' @param bed_file BED file with exon/intron junction annotations
#' @param output_file Path to output result file
#' @param min_reads Minimum number of supporting reads to report a novel junction (default 5)
#' @param return_granges Logical; if TRUE, return result as a GRanges object with supporting reads in list-columns
#' @return If return_granges is TRUE, a GRanges object is returned; otherwise nothing is returned
#' @export

mySplice <- function(maf_file, bam_list, bed_file, output_file, min_reads = 5, return_granges = TRUE) {
  # Load Python function using reticulate
  py_script <- system.file("python", "in_silico_ns.py", package = "mySpliceR")
  reticulate::source_python(py_script)

  # Call the Python function
  result_list = detect_novel_junctions(
    maf_file = maf_file,
    bam_list = bam_list,
    bed_file = bed_file,
    output_file = output_file,
    min_reads = min_reads
  )

# If GRanges output requested, convert Python list to GRanges
  if (return_granges) {
    if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
      stop("Package 'GenomicRanges' is required to return GRanges object")
    }

    variant_keys <- names(result_list)
    variant_data <- result_list[variant_keys]

    gr <- GenomicRanges::GRanges(
      seqnames = sapply(variant_data, function(x) x$chrom),
      ranges = IRanges::IRanges(
        start = sapply(variant_data, function(x) x$pos),
        end   = sapply(variant_data, function(x) x$pos)
      ),
      strand = "*",
      sample = sapply(variant_data, function(x) x$sample),
      ref = sapply(variant_data, function(x) x$ref),
      alt = sapply(variant_data, function(x) x$alt),
      num_supporting_reads = sapply(variant_data, function(x) x$num_supporting_reads),
      supporting_reads = S4Vectors::SimpleList(lapply(variant_data, function(x) {
        data.frame(
          read_name = sapply(x$supporting_reads, function(r) r$read_name),
          sequence = sapply(x$supporting_reads, function(r) r$sequence),
          junction_start = sapply(x$supporting_reads, function(r) r$junction_start),
          junction_end = sapply(x$supporting_reads, function(r) r$junction_end),
          stringsAsFactors = FALSE
        )
      }))
    )
    return(gr)
  } else {
    message("Mis-splice detection completed. Output written to: ", output_file)
    invisible(NULL)
  }
}

