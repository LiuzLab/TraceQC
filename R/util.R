
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

#' Get absolute path of a file.
#'
#' @param f A relative or absolute file path.
#'
#' @return It returns an absolute path for a file.
get_abspath <- function(f) {
  if(is.null(f)) {
    stop("NULL shouldn't be an argument.")
  }
  file.path(normalizePath(dirname(f)), basename(f))
}


#' Split a string by a fixed length and joined with <br> HTML tag.
#'
#' @param s The input string
#' @param len The fixed length
#'
#' @import stringr
#'
#' @return A string splited by the `len' then joined by `<br/>'
#'
#' @export
seq_split <- function(s, len = 50) {
  st <- seq(1, str_length(s), by=len)
  ed <- st+len
  str_replace_all(paste0(str_sub(s, st, ed-1), collapse = "<br/>"), "-", "_")
}


#' Installation necessary external software for TraceQC
#'
#' @export
#'
#' @importFrom fastqcr fastqc_install
#' @importFrom reticulate py_install
#'
#' @examples
#' \dontrun{
#' library(TraceQC)
#' install_external_packages()
#' }
install_external_packages <- function() {
  fastqc_install()
  py_install(packages = c("pandas", "biopython", "tqdm"))
  message("All external packages has been installed.")
}

