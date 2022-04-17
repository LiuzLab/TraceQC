#' Function for a sequence alignment between the reference file and sample.
#'
#' The function is an wrapper of a python function which performs
#' a global pairwise sequence alignment by biopython package.
#'
#' @param input_file A FASTQ file path (required).
#' @param ref_file A path of a reference sequence file (required).
#' @param output_file The output path. An output of the alignment will be
#' stored at the path (default: "aligned_reads.txt").
#' @param python_path The path to Python. (default: "python3").
#' @param match The score for a correct basepair matching (default: 2).
#' @param mismatch The penalty score for a basepair mismatching (default: -2).
#' @param gapopen The gap opening score for the alignment (default: -6).
#' @param gapextension The gap extension score for the alignment (default: -0.1).
#' @param ncores The number of cores for the parallel processing.
#' @param penalize_end_gaps Whether penalize the end gap, the value should take 0/1 (default: 1).
#' @param return_df A logical argument to report what type of output will
#' be created from the function.
#'
#' @importFrom readr read_tsv
#' @return It returns a data frame of the alignment result
#' if `return_df' is `T' and `NULL' otherwise.
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
                               output_file="aligned_reads.txt",
                               python_path = "python3",
                               match=2,
                               mismatch=-2,
                               gapopen=-6,
                               gapextension=-0.1,
                               ncores = 4,
                               penalize_end_gaps=1,
                               return_df = FALSE) {
  path <- system.file("py", "alignment.py", package="TraceQC")
  system(sprintf("%s %s --input %s --ncores %s --ref %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s",
                 python_path,path,input_file,ncores,ref_file,output_file,match,mismatch,gapopen,gapextension))}

#' Function for a sequence alignment between the reference file and sample for 10x data.
#'
#' The function is an wrapper of a python function which performs
#' a global pairwise sequence alignment by biopython package.
#'
#' @param input_file A FASTQ file path (required).
#' @param ref_file A path of a reference sequence file (required).
#' @param output_file The output path. An output of the alignment will be
#' stored at the path (default: "aligned_reads.txt").
#' @param python_path The path to Python. (default: "python3").
#' @param match The score for a correct basepair matching (default: 2).
#' @param mismatch The penalty score for a basepair mismatching (default: -2).
#' @param gapopen The gap opening score for the alignment (default: -6).
#' @param gapextension The gap extension score for the alignment (default: -0.1).
#' @param ncores The number of cores for the parallel processing.
#' @param penalize_end_gaps Whether penalize the end gap, the value should take 0/1 (default: 1).
#'
#' @importFrom readr read_tsv
#' @return It returns a data frame of the alignment result
#' if `return_df' is `T' and `NULL' otherwise.
#' @export
#'
#' @examples
#' input_file <- system.file("extdata/test_data/raw_sequence","carlin_example.bam",package="TraceQC")
#' output_file <- "./aligned_reads.txt"
#' sequence_alignment_for_10x(input_file=input_file,ref_file=ref_file,
#'                            output_file=output_file)

sequence_alignment_for_10x <- function(input_file,
                                       ref_file,
                                       output_file="aligned_reads.txt",
                                       python_path = "python3",
                                       match=2,
                                       mismatch=-2,
                                       gapopen=-6,
                                       gapextension=-0.1,
                                       penalize_end_gaps=1,
                                       ncores = 4) {
  path <- system.file("py", "alignment_for_10x.py", package="TraceQC")
  system(sprintf("%s %s --input %s --ncores %s --ref %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s --penalize_end_gaps %s",
                 python_path,path,input_file,ncores,ref_file,output_file,match,mismatch,gapopen,gapextension,penalize_end_gaps))}

#' Function for finding threshold of sequence alignment. The function randomly
#' permutate certain percentage reference sequence and perform global alignment
#' with the original reference sequence. By use the permutated sequence alignment
#' score, users can filter the TraceQC alignment result.
#'
#' @param ref_file A path of a reference sequence file (required).
#' @param output_file The output path. An output of the alignment will be
#' stored at the path (default: "alignment_threshold.txt").
#' @param python_path The path to Python. (default: "python3").
#' @param match The score for a correct basepair matching (default: 2).
#' @param mismatch The penalty score for a basepair mismatching (default: -2).
#' @param gapopen The gap opening score for the alignment (default: -6).
#' @param gapextension The gap extension score for the alignment (default: -0.1).
#' @param penalize_end_gaps Whether penalize the end gap, the value should take 0/1 (default: 1).
#'
#' @importFrom readr read_tsv
#' @return It returns a data frame of the alignment result with randomly permutated sequence.
#' @export
#' @examples
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' sequence_permutation(ref_file,out="./seq_permutate.txt")

sequence_permutation <- function(ref_file,
                                python_path = "python3",
                                match=2,
                                mismatch=-2,
                                gapopen=-6,
                                gapextension=-0.1,
                                penalize_end_gaps=1,
                                read_length=0,
                                permutate_percent=seq(0,1,length.out=101),
                                n=2,
                                output_file="alignment_threshold.txt") {
  permutate_file <- tempfile()
  write(paste(permutate_percent,collapse=" "), file=permutate_file)
  path <- system.file("py", "sequence_alignment_threshold.py", package="TraceQC")
  system(sprintf("%s %s --ref %s --n %s --permutate_percent %s --read_length %s --output %s --match %s --mismatch %s --gapopen %s --gapextension %s --penalize_end_gaps %s",
                 python_path,path,ref_file,n,permutate_file,read_length,output_file,match,mismatch,gapopen,gapextension,penalize_end_gaps))
  file.remove(permutate_file)}


#' Parsing reference sequence file
#' @param ref_file A path of a reference sequence file.
#' @return A list with those four elements.
#' \itemize{
#'   \item `refseq': The reference sequence.
#'   \item `regions': Detailed information about the reference sequence.
#' }
#' @export
#' @examples
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' ref <- parse_ref_file(ref_file)
#' 
parse_ref_file <- function(ref_file) {
  ref <- readLines(ref_file)
  refseq <- ref[1]
  regions <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions) %>%
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(.data$start),
           end=strtoi(.data$end)) %>%
    mutate(region=as.character(.data$region))
  list(refseq=refseq,regions=regions)
}
