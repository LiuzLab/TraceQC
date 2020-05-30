#' Function to locate where the Quality Control file is stored.
#'
#' @param input_file A path of FASTQ file.
#' @param qc_dir The directory path where Quality Control files are.
#'
#' @return A path to the corresponded zipped Quality Control file.
#' @export
get_qcpath <- function(input_file, qc_dir) {
  qc_path <- sub(".fastq", "_fastqc.zip", basename(input_file))
  file.path(qc_dir, qc_path)
}
