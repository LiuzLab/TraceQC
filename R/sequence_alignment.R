#' Function for a sequence alignment between the reference file and sample.
#'
#' This function serves as a wrapper for a Python script that performs global pairwise
#' sequence alignment using the Biopython package. It aligns sequencing reads from a FASTQ file
#' against a reference sequence and outputs the results to a specified file.
#'
#' @param input_file A character string specifying the path to the input FASTQ file (required).
#' @param ref_file A character string specifying the path to the reference sequence file (required).
#' @param output_file A character string specifying the output file path (default: "aligned_reads.txt").
#' @param python_path A character string specifying the path to the Python executable (default: "python3").
#' @param match Numeric value for the score of a correct basepair match (default: 2).
#' @param mismatch Numeric value for the penalty of a basepair mismatch (default: -2).
#' @param gapopen Numeric value for the gap opening penalty (default: -6).
#' @param gapextension Numeric value for the gap extension penalty (default: -0.1).
#' @param ncores Integer specifying the number of CPU cores for parallel processing (default: 4).
#' @param penalize_end_gaps Integer (0 or 1) indicating whether to penalize end gaps (default: 1).
#' @param return_df Logical indicating whether to return the alignment results as a data frame (default: FALSE).
#'
#' @importFrom readr read_tsv
#' @return If `return_df` is TRUE, returns a data frame containing the alignment results.
#'   Otherwise, returns NULL invisibly. The alignment results are always written to `output_file`.
#' @export
#'
#' @examples
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' ref <- parse_ref_file(ref_file)
#' input_file <- system.file("extdata/test_data/raw_sequence","hgRNA_example.fastq.gz",package="TraceQC")
#' output_file <- "./aligned_reads.txt"
#' sequence_alignment(input_file=input_file,ref_file=ref_file,
#'                    output_file=output_file)
#'
sequence_alignment <- function(input_file,
                               ref_file,
                               output_file = "aligned_reads.txt",
                               python_path = "python3",
                               match = 2,
                               mismatch = -2,
                               gapopen = -6,
                               gapextension = -0.1,
                               ncores = 4,
                               penalize_end_gaps = 1,
                               return_df = FALSE) {
  # Input validation
  if (!file.exists(input_file)) stop("Input file does not exist: ", input_file)
  if (!file.exists(ref_file)) stop("Reference file does not exist: ", ref_file)
  if (ncores < 1) stop("ncores must be a positive integer")

  # Locate the Python script within the package
  py_script <- system.file("py", "alignment.py", package = "TraceQC")
  if (py_script == "") stop("Python alignment script not found in package")

  # Construct the command to run the Python script
  cmd <- sprintf(
    "%s %s --input %s --ncores %s --ref %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s --penalize_end_gaps %s",
    python_path, py_script, shQuote(input_file), ncores, shQuote(ref_file),
    shQuote(output_file), match, mismatch, gapopen, gapextension, penalize_end_gaps
  )

  # Execute the Python script
  exit_code <- system(cmd)
  if (exit_code != 0) stop("Python alignment script failed with exit code: ", exit_code)

  # Optionally return the results as a data frame
  if (return_df) {
    if (!file.exists(output_file)) stop("Output file was not created: ", output_file)
    return(read_tsv(output_file, show_col_types = FALSE))
  }
  invisible(NULL)
}

#' Function for a sequence alignment between the reference file and sample for 10x data.
#'
#' This function is a wrapper for a Python script that performs global pairwise sequence alignment
#' using the Biopython package, specifically tailored for 10x Genomics single-cell data (BAM format).
#' It aligns sequencing reads against a reference sequence and outputs the results.
#'
#' @param input_file A character string specifying the path to the input BAM file (required).
#' @param ref_file A character string specifying the path to the reference sequence file (required).
#' @param output_file A character string specifying the output file path (default: "aligned_reads.txt").
#' @param python_path A character string specifying the path to the Python executable (default: "python3").
#' @param match Numeric value for the score of a correct basepair match (default: 2).
#' @param mismatch Numeric value for the penalty of a basepair mismatch (default: -2).
#' @param gapopen Numeric value for the gap opening penalty (default: -6).
#' @param gapextension Numeric value for the gap extension penalty (default: -0.1).
#' @param penalize_end_gaps Integer (0 or 1) indicating whether to penalize end gaps (default: 1).
#' @param ncores Integer specifying the number of CPU cores for parallel processing (default: 4).
#'
#' @importFrom readr read_tsv
#' @return Always returns NULL invisibly. Results are written to `output_file`.
#' @export
#'
#' @examples
#' input_file <- system.file("extdata/test_data/raw_sequence","carlin_example.bam",package="TraceQC")
#' output_file <- "./aligned_reads.txt"
#' sequence_alignment_for_10x(input_file=input_file,ref_file=ref_file,
#'                            output_file=output_file)

sequence_alignment_for_10x <- function(input_file,
                                       ref_file,
                                       output_file = "aligned_reads.txt",
                                       python_path = "python3",
                                       match = 2,
                                       mismatch = -2,
                                       gapopen = -6,
                                       gapextension = -0.1,
                                       penalize_end_gaps = 1,
                                       ncores = 4) {
  # Input validation
  if (!file.exists(input_file)) stop("Input file does not exist: ", input_file)
  if (!file.exists(ref_file)) stop("Reference file does not exist: ", ref_file)
  if (ncores < 1) stop("ncores must be a positive integer")

  # Locate the Python script within the package
  py_script <- system.file("py", "alignment_for_10x.py", package = "TraceQC")
  if (py_script == "") stop("Python alignment script for 10x not found in package")

  # Construct the command to run the Python script
  cmd <- sprintf(
    "%s %s --input %s --ncores %s --ref %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s --penalize_end_gaps %s",
    python_path, py_script, shQuote(input_file), ncores, shQuote(ref_file),
    shQuote(output_file), match, mismatch, gapopen, gapextension, penalize_end_gaps
  )

  # Execute the Python script
  exit_code <- system(cmd)
  if (exit_code != 0) stop("Python alignment script failed with exit code: ", exit_code)

  invisible(NULL)
}#' Function for finding threshold of sequence alignment via permutation testing.
#'
#' This function performs permutation testing to determine optimal alignment score thresholds.
#' It randomly permutes portions of the reference sequence and aligns them back to the original,
#' generating a distribution of alignment scores that can be used to filter true alignments
#' from noise in downstream analysis.
#'
#' @param ref_file A character string specifying the path to the reference sequence file (required).
#' @param python_path A character string specifying the path to the Python executable (default: "python3").
#' @param match Numeric value for the score of a correct basepair match (default: 2).
#' @param mismatch Numeric value for the penalty of a basepair mismatch (default: -2).
#' @param gapopen Numeric value for the gap opening penalty (default: -6).
#' @param gapextension Numeric value for the gap extension penalty (default: -0.1).
#' @param penalize_end_gaps Integer (0 or 1) indicating whether to penalize end gaps (default: 1).
#' @param read_length Integer specifying the read length for simulation (default: 0, uses full ref length).
#' @param permutate_percent Numeric vector of permutation percentages to test (default: seq(0,1,length.out=101)).
#' @param n Integer specifying the number of permutations per percentage (default: 2).
#' @param output_file A character string specifying the output file path (default: "alignment_threshold.txt").
#'
#' @importFrom readr read_tsv
#' @return Always returns NULL invisibly. Results are written to `output_file`.
#' @export
#' @examples
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' sequence_permutation(ref_file, output_file = "./seq_permutate.txt")

sequence_permutation <- function(ref_file,
                                 python_path = "python3",
                                 match = 2,
                                 mismatch = -2,
                                 gapopen = -6,
                                 gapextension = -0.1,
                                 penalize_end_gaps = 1,
                                 read_length = 0,
                                 permutate_percent = seq(0, 1, length.out = 101),
                                 n = 2,
                                 output_file = "alignment_threshold.txt") {
  # Input validation
  if (!file.exists(ref_file)) stop("Reference file does not exist: ", ref_file)
  if (n < 1) stop("n must be a positive integer")
  if (any(permutate_percent < 0 | permutate_percent > 1)) stop("permutate_percent must be between 0 and 1")

  # Create a temporary file for permutation percentages
  perm_file <- tempfile()
  on.exit(file.remove(perm_file), add = TRUE)
  writeLines(paste(permutate_percent, collapse = " "), perm_file)

  # Locate the Python script within the package
  py_script <- system.file("py", "sequence_alignment_threshold.py", package = "TraceQC")
  if (py_script == "") stop("Python permutation script not found in package")

  # Construct the command to run the Python script
  cmd <- sprintf(
    "%s %s --ref %s --n %s --permutate_percent %s --read_length %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s --penalize_end_gaps %s",
    python_path, py_script, shQuote(ref_file), n, shQuote(perm_file), read_length,
    shQuote(output_file), match, mismatch, gapopen, gapextension, penalize_end_gaps
  )

  # Execute the Python script
  exit_code <- system(cmd)
  if (exit_code != 0) stop("Python permutation script failed with exit code: ", exit_code)

  invisible(NULL)
}


#' Parsing reference sequence file
#'
#' This function reads and parses a reference sequence file used in TraceQC.
#' The file format consists of the full construct sequence on the first line,
#' followed by lines defining regions (e.g., "target 24 140").
#'
#' @param ref_file A character string specifying the path to the reference sequence file.
#' @return A list containing:
#'   \itemize{
#'     \item `refseq`: The full reference sequence as a character string.
#'     \item `regions`: A data frame with columns `region`, `start`, and `end` defining sequence regions.
#'   }
#' @export
#' @examples
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' ref <- parse_ref_file(ref_file)
#'
parse_ref_file <- function(ref_file) {
  # Input validation
  if (!file.exists(ref_file)) stop("Reference file does not exist: ", ref_file)

  # Read all lines from the reference file
  ref_lines <- readLines(ref_file, warn = FALSE)

  # Extract the reference sequence (first line)
  refseq <- ref_lines[1]

  # Parse region definitions (subsequent lines)
  if (length(ref_lines) > 1) {
    regions <- strsplit(ref_lines[2:length(ref_lines)], split = " ")
    regions <- do.call(rbind, regions) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      setNames(c("region", "start", "end")) %>%
      mutate(
        start = as.integer(.data$start),
        end = as.integer(.data$end),
        region = as.character(.data$region)
      )
  } else {
    # No regions defined, create empty data frame
    regions <- data.frame(region = character(), start = integer(), end = integer())
  }

  # Return as a list
  list(refseq = refseq, regions = regions)
}
