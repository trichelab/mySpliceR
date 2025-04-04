#' Detect mis-splicing events based on somatic mutations
#'
#' @param maf_file Path to the input MAF file
#' @param bam_list Path to the TSV file listing sample ID and BAM paths
#' @param bed_file BED file with exon/intron junction annotations
#' @param output_file Path to output result file
#' @param min_reads Minimum number of supporting reads to report a novel junction (default 5)
#' @export
mySplice <- function(maf_file, bam_list, bed_file, output_file, min_reads = 5) {
  # Load Python function using reticulate
  py_script <- system.file("python", "in_silico_ns.py", package = "mySpliceR")
  reticulate::source_python(py_script)

  # Call the Python function
  detect_novel_junctions(
    maf_file = maf_file,
    bam_list = bam_list,
    bed_file = bed_file,
    output_file = output_file,
    min_reads = min_reads
  )

  message("Mis-splice detection completed. Output written to: ", output_file)
}

