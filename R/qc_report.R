

#' Create a HTML QC report
#'
#' @param input_file A path of a FASTQ file.
#' @param ref_file A path of a text file which contains construct information.
#' @param fastqc_dir A directory path to store FastQC files. It will be omiited
#' if it is set as NULL
#' @param output_path A file path to save the HTML report. It will be omiited
#' if it is set as NULL
#' @param ncores The number of cores for parallel processing.
#' @param title The title of the report.
#' @param preview It will open a preview to a web-browser if the argument is set to `TRUE'
#'
#' @return A TraceQC object for the input will be returned.
#'
#' @imporFrom rmarkdown render
#'
#' @export
#'
#' @examples
#' library(TraceQC)
#' obj <- generate_qc_report(
#'   input_file = system.file("extdata", "test_data",
#'                            "fastq", "example_14d.fastq", package="TraceQC"),
#'   ref_file = system.file("extdata", "test_data",
#'                          "ref", "ref.txt", package="TraceQC"),
#'   preview = FALSE,
#'   title = "TraceQC report for example_14d.fastq",
#'   ncores=1
#'   )
#' summary(obj)
#'
generate_qc_report <-
  function(
    input_file = NULL,
    ref_file = NULL,
    fastqc_dir = NULL,
    output_path = NULL,
    ncores = 4,
    title = "TraceQC report",
    preview = FALSE
  ) {

    if(is.null(input_file)) {
      stop("input_file shouldn't be NULL.")
    }

    if(!file.exists(input_file)) {
      stop("input_file is't existed.")
    }

    if(is.null(ref_file)) {
      stop("ref_file shouldn't be NULL.")
    }

    if(!file.exists(ref_file)) {
      stop("ref_file is't existed.")
    }

    if(is.null(fastqc_dir)) {
      fastqc_dir <- tempdir()
    }

    if(is.null(output_path)) {
      fastq <- tempfile(fileext = ".html")
    }
    template_path <- system.file(
      "Rmd",
      "TraceQC-template.Rmd",
      package = "TraceQC"
    )

    rds_path <- tempfile(fileext = ".rds")
    knitr_params <- list()
    knitr_params$debug <- FALSE
    knitr_params$input_file <- input_file
    knitr_params$ref_file <- ref_file
    knitr_params$fastqc_dir <- fastqc_dir
    knitr_params$ncores <- ncores
    knitr_params$date <- Sys.Date()
    knitr_params$rds_path <- rds_path
    knitr_params$set_title <- title

    rmdout_path <- rmarkdown::render(
      input = template_path,
      output_format = "html_document",
      output_file = output_path,
      params = knitr_params
    )

    if(preview)
      browseURL(paste0("file://", rmdout_path))

    readRDS(rds_path)
  }


