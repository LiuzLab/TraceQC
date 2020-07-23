
#' Function to locate where the Quality Control file is stored.
#'
#' @param input_file A path of FASTQ file.
#' @param qc_dir The directory path where Quality Control files are.
#'
#' @return A path to the corresponded zipped Quality Control file.
#' @export
get_qcpath <- function(input_file, qc_dir) {
  qc_path <- strsplit(basename(input_file), split = "\\.")[[1]][1]
  qc_path <-  paste0(qc_path, "_fastqc.zip")
  path <- file.path(qc_dir, qc_path)
  if(!file.exists(path)) {
    stop(
    paste("A FastQC file for the input file was not found.",
          "Please be aware that TraceQC is not guranteed to run for FASTQ files that contain period before the file extension.",
          "For example, Running a sample with sample.02.fastq.gz may occur some errors.",
          sep = "\n")
    )
  }
  path
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


#' Running FASTQC for a file.
#'
#' @param input_file A FASTQ file path
#' @param threads The number of threads will be used in FASTQC
#'
#' @return It returns a path to the QC file.
#' @export
#'
#' @examples
#' input_file <- system.file("extdata", "test_data",
#'                           "fastq", "example.fastq.gz", package="TraceQC")
#' fastqc.file(input_file)
fastqc.file <- function(input_file,
                        threads = 4) {
  tmp_dir <- tempdir()
  fq_dir <- file.path(tmp_dir, basename(input_file), "fq")
  qc_dir <- file.path(tmp_dir, basename(input_file), "qc")
  dir.create(fq_dir, recursive=T)
  dir.create(qc_dir, recursive=T)

  fq_path <- file.path(fq_dir, basename(input_file))
  if(!file.exists(fq_path))
    file.symlink(input_file, fq_path)
  ref_file <- system.file("extdata", "test_data", "ref",
                          "ref.txt", package="TraceQC")

  fastqcr::fastqc(fq.dir = fq_dir, qc.dir=qc_dir, threads = threads)
  qc_dir
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
