
#' Installation necessary external software for TraceQC
#'
#' @return
#' @export
#'
#' @importFrom fastqcr fastqc_install
#' @importFrom reticulate py_install
#'
#' @examples
#' library(TraceQC)
#' install_external_packages()
install_external_packages <- function() {
  fastqc_install()
  py_install(packages = c("pandas", "biopython"))
  message("All external packages has been installed.")
}

