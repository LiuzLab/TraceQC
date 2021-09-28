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

