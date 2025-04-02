test_that("split_maf produces the correct number of files", {
  # Setup test input
  input_file <- system.file("extdata", "example.maf", package = "mySpliceR")
  output_dir <- tempfile("split_test_")
  max_files <- 3

  # Run split
  split_maf(input_file, output_dir, max_files)

  # Check output
  out_files <- list.files(output_dir, pattern = "\\.maf\\.\\d+$", full.names = TRUE)
  expect_true(length(out_files) == max_files)

  # Count rows in the original file (excluding header)
  original_rows <- length(readLines(input_file)) - 1

  # Check each file has the header
  header_line <- readLines(input_file, n = 1)
  for (f in out_files) {
    lines <- readLines(f)
    expect_equal(lines[1], header_line)
  }
  # Sum up all data rows in the split files (excluding header)
  split_rows <- sum(sapply(out_files, function(f) {
    length(readLines(f)) - 1
  }))

  expect_equal(split_rows, original_rows)

})

