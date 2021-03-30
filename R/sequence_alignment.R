#' Function for a sequence alignment between the reference file and sample.
#'
#' The function is an wrapper of a python function which performs
#' a global pairwise sequence alignment by biopython package.
#'
#' @param input_file A FASTQ file path
#' @param ref_file A path of a reference sequence file.
#' @param output_file The output path. An output of the alignment will be
#' stored at the path.
#' @param match The score for a correct basepair matching.
#' @param mismatch The penalty score for a basepair mismatching.
#' @param gapopen The gap opening score for the alignment.
#' @param gapextension The gap extension score for the alignment.
#' @param return_df A logical argument to report what type of output will
#' be created from the function.
#' @param ncores The number of cores for the parallel processing
#'
#' @importFrom reticulate source_python
#' @importFrom readr read_tsv
#' @importFrom tictoc tic toc
#' @return It returns a data frame of the alignment result
#' if `return_df' is `T' and `NULL' otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' library(TraceQC)
#' input_file <- system.file("extdata", "test_data",
#'                           "fastq", "example_small.fastq.gz", package="TraceQC")
#' ref_file <- system.file("extdata", "test_data", "ref",
#'                         "ref.txt", package="TraceQC")
#' output_file <- tempfile()
#' sequence_alignment(input_file=input_file,
#'                    ref_file=ref_file,
#'                    output_file=output_file,
#'                    return_df=TRUE)
#' }
sequence_alignment <- function(input_file,
                               ref_file,
                               output_file="aligned_reads.txt",
                               match=2,
                               mismatch=-2,
                               gapopen=-6,
                               gapextension=-0.1,
                               ncores = 4,
                               return_df = FALSE) {
  args <- list("input"=get_abspath(input_file),
               "reference"=get_abspath(ref_file),
               "output"=get_abspath(output_file),
               "match"=match,"mismatch"=mismatch,"gapopen"=gapopen,
               "gapextension"=gapextension,
               "ncores"=ncores)
  source_python(system.file("py", "alignment.py", package="TraceQC"))
  print(ls(parent.frame()))
  module <- reticulate::import_from_path("alignment", system.file("py", package = "TraceQC"))


  tic("Alignment")
  message(paste0("Running an alignment between ", ref_file,
                 " and ", input_file, "."))
  module$alignment(args)
  toc()

  if(return_df) {
    read_tsv(output_file)
  } else {
    return(NULL)
  }
}

#' Function for a sequence alignment between the reference file and sample for 10x data.
#'
#' The function is an wrapper of a python function which performs
#' a global pairwise sequence alignment by biopython package.
#'
#' @param input_file  A file path of possorted_genome_bam.bam out put by cellranger
#' @param ref_file A path of a reference sequence file.
#' @param output_file The output path. An output of the alignment will be
#' stored at the path.
#' @param match The score for a correct basepair matching.
#' @param mismatch The penalty score for a basepair mismatching.
#' @param gapopen The gap opening score for the alignment.
#' @param gapextension The gap extension score for the alignment.
#' @param return_df A logical argument to report what type of output will
#' be created from the function.
#' @param ncores The number of cores for the parallel processing
#'
#' @importFrom reticulate source_python
#' @importFrom readr read_tsv
#' @importFrom tictoc tic toc
#' @return It returns a data frame of the alignment result
#' if `return_df' is `T' and `NULL' otherwise.
#' @export
#'
#' @examples
sequence_alignment_for_10x <- function(input_file,
                               ref_file,
                               output_file="aligned_reads.txt",
                               match=2,
                               mismatch=-2,
                               gapopen=-6,
                               gapextension=-0.1,
                               penalize_end_gaps=0,
                               ncores = 4,
                               return_df = FALSE) {
  args <- list("input"=get_abspath(input_file),
               "reference"=get_abspath(ref_file),
               "output"=get_abspath(output_file),
               "match"=match,"mismatch"=mismatch,"gapopen"=gapopen,
               "gapextension"=gapextension,"penalize_end_gaps"=penalize_end_gaps,
               "ncores"=ncores)
  source_python(system.file("py", "alignment_for_10x.py", package="TraceQC"))
  print(ls(parent.frame()))
  module <- reticulate::import_from_path("alignment_for_10x", system.file("py", package = "TraceQC"))


  tic("Alignment")
  message(paste0("Running an alignment between ", ref_file,
                 " and ", input_file, "."))
  module$alignment(args)
  toc()

  if(return_df) {
    read_tsv(output_file)
  } else {
    return(NULL)
  }
}

#' Function for finding threshold of sequence alignment. The function randomly
#' permutate certain percentage reference sequence and perform global alignment
#' with the original reference sequence. By use the permutated sequence alignment
#' score, users can filter the TraceQC alignment result.
#'
#' @param ref_file A path of a reference sequence file.
#' @param output_file The output path. An output dataframe will be
#' stored at the path.
#' @param match The score for a correct basepair matching.
#' @param mismatch The penalty score for a basepair mismatching.
#' @param gapopen The gap opening score for the alignment.
#' @param gapextension The gap extension score for the alignment.
#' @param n number of random permutation used for each percentage
#' @param corrupted percentage The number of cores for the parallel processing
#'
#' @importFrom reticulate source_python
#' @importFrom readr read_tsv
#' @return It returns a data frame of the alignment result
#' @export

sequnce_alignment_threshold <- function(ref_file,
                                        match=2,
                                        mismatch=-2,
                                        gapopen=-6,
                                        gapextension=-0.1,
                                        penalize_end_gaps=1,
                                        read_length=0,
                                        permutate_percent=c(0.1,0.2,0.3,0.4,0.5,0.6),
                                        n=100,
                                        output_file="alignment_threshold.txt") {

    args <- list("reference"=get_abspath(ref_file),
                 "output_file"=get_abspath(output_file),
                 "match"=match,"mismatch"=mismatch,"gapopen"=gapopen,
                 "penalize_end_gaps"=penalize_end_gaps,"read_length"=read_length,
                 "gapextension"=gapextension,"n"=n,
                 "permutate_percent"=permutate_percent)
    # source_python(system.file("py", "sequence_alignment_threshold.py", package="TraceQC"))
    module <- reticulate::import_from_path("sequence_alignment_threshold", system.file("py", package = "TraceQC"))
    module$alignment_score_threshold(args)
    read_tsv(output_file)
}
