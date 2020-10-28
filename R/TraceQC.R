#' Creating an input object for TraceQC
#'
#' @param input_file A FASTQ file path.
#' @param ref_file A path of a reference sequence file.
#' @param fastqc_file A path of a FASTQC file.
#' @param alignment_output_file A path to store alignment output file.
#' This will be `NULL' if the output is temporally stored.
#' @param ncores The number of cores for the parallel processing.
#'
#' @return It will return a list from `create_input_object_with_alignment'
#' @export
#'
#' @examples
#' library(TraceQC)
#' library(fastqcr)
#' input_file <- system.file("extdata", "test_data",
#'                           "fastq", "example_small.fastq.gz", package="TraceQC")
#' ref_file <- system.file("extdata", "test_data", "ref",
#'                         "ref.txt", package="TraceQC")
#' qc_dir <- fastqc.file(input_file)
#'
#' input_qc_path <- get_qcpath(input_file, qc_dir)
#'
#' obj <- TraceQC(input_file = input_file,
#'                ref_file = ref_file,
#'                fastqc_file = input_qc_path,
#'                ncores=1)
#'
#' obj$refseq
TraceQC <-
  function(input_file,
           ref_file,
           fastqc_file=NULL,
           alignment_output_file=NULL,
           ncores=4) {

  if(is.null(alignment_output_file)) {
    alignment_output_file <- tempfile()
  }

  sequence_alignment(
    input_file = input_file,
    ref_file = ref_file,
    output_file = alignment_output_file,
    ncores = ncores
  )
  create_TraceQC_object(
    alignment_output_file,
    ref_file,
    fastqc_file,
    ncores
    )
}


#' Creating an input object for TraceQC with an alignment result.
#'
#' @param aligned_reads_file A path to store alignment output file.
#' @param ref_file A path of a reference sequence file.
#' @param fastqc_file A path of a FASTQC file.
#' @param ncores The number of cores for the parallel processing.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @importFrom readr read_tsv
#' @importFrom fastqcr qc_read
#' @importFrom stats setNames
#'
#' @return A list with those four elements.
#' \itemize{
#'   \item `aligned_reads': an data frame that contains alignment information.
#'   \item `refseq': The reference sequence.
#'   \item `regions': Detailed information about the reference sequence.
#'   \item `qc': A list of tibbles containing the QC data.
#' }
#' @export
#'
create_TraceQC_object <-
  function(aligned_reads_file,
           ref_file,
           fastqc_file,
           use_CPM = TRUE,
           alignment_score_threshold=0,
           abundance_threshold=0,
           ncores=1) {

    aligned_reads <- read_tsv(aligned_reads_file)

    ref <- readLines(ref_file)
    refseq <- ref[1]
    regions <- strsplit(ref[2:length(ref)],split=" ")
    regions <- do.call(rbind,regions) %>%
      as.data.frame() %>%
      setNames(c("region","start","end")) %>%
      mutate(start=strtoi(.data$start),
             end=strtoi(.data$end)) %>%
      mutate(region=as.character(.data$region))

    qc <- NULL
    if(!is.null(fastqc_file)) {
      qc <- qc_read(fastqc_file)
    }

    obj <- list(aligned_reads=aligned_reads,
                refseq=refseq,
                regions=regions,
                qc=qc)

    message("Running mutation event identification.")
    tic("mutation event identification")
    obj$mutation <- seq_to_character(obj,
                                     ncores,
                                     use_CPM,
                                     alignment_score_threshold,
                                     abundance_threshold)
    toc()
    obj
  }
