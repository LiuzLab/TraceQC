#' Get absolute path of a file.
#'
#' This utility function normalizes a file path to its absolute form, expanding
#' relative paths and resolving symbolic links. It's useful for ensuring consistent
#' file references across different working directories.
#'
#' @param f A character string representing a file path (relative or absolute).
#'
#' @return A character string containing the absolute path to the file.
#' @export
get_abspath <- function(f) {
  if (is.null(f)) {
    stop("Input path cannot be NULL")
  }
  if (!is.character(f) || length(f) != 1) {
    stop("Input must be a single character string")
  }
  file.path(normalizePath(dirname(f)), basename(f))
}


#' Split a string by a fixed length and join with HTML break tags.
#'
#' This function splits a long string (typically a DNA sequence) into chunks of
#' specified length and joins them with HTML <br/> tags for display purposes,
#' such as in R Markdown documents or Shiny apps. It also converts dashes to
#' underscores for better HTML rendering.
#'
#' @param s A character string to split (e.g., a DNA sequence).
#' @param len An integer specifying the chunk length (default: 50).
#'
#' @import stringr
#' @return A character string with chunks separated by <br/> tags.
#' @export
seq_split <- function(s, len = 50) {
  if (!is.character(s) || length(s) != 1) {
    stop("Input must be a single character string")
  }
  if (!is.numeric(len) || len <= 0) {
    stop("len must be a positive number")
  }

  # Calculate chunk boundaries
  st <- seq(1, str_length(s), by = len)
  ed <- pmin(st + len - 1, str_length(s))

  # Extract chunks and join with HTML breaks
  chunks <- str_sub(s, st, ed)
  result <- paste0(chunks, collapse = "<br/>")

  # Replace dashes with underscores for HTML compatibility
  str_replace_all(result, "-", "_")
}

