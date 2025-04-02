## R/mySplice.R
mySplice <- function(maf_file, bam_list, bed_file, output_file, min_reads = 5) {
  # Load the Python module using reticulate
  py_script <- system.file("python", "in_silico_ns.py", package = "MySpliceR")
  reticulate::source_python(py_script)

  # Call Python function (assumes it exposes `detect_novel_junctions`)
  detect_novel_junctions(maf_file, bam_list, bed_file, output_file, min_reads)
  
  message("Processing complete. Output written to:", output_file)
}
