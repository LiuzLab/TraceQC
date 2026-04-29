#' Get read count per UMI.
#'
#' This function calculates the number of reads associated with each UMI (Unique Molecular Identifier)
#' within each cell barcode. This is used for quality control and filtering in single-cell data.
#'
#' @param df A data frame of aligned reads containing 'UB' (UMI) and 'CB' (cell barcode) columns.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A data frame with columns: UB, CB, read_count_per_UMI.
#' @export
#'
get_read_count_per_UB <- function(df) {
  # Input validation
  if (!all(c("UB", "CB") %in% names(df))) {
    stop("Data frame must contain 'UB' and 'CB' columns")
  }

  df %>%
    group_by(UB, CB) %>%
    summarise(read_count_per_UMI = n(), .groups = "drop")
}


#' Get UMI count per cell.
#'
#' This function calculates the number of unique UMIs per cell barcode, which is a measure
#' of the sequencing depth and complexity for each cell.
#'
#' @param df A data frame of aligned reads containing 'UB' (UMI) and 'CB' (cell barcode) columns.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A data frame with columns: CB, UMI_per_CB.
#' @export
#'
get_UMI_count_per_CB <- function(df) {
  # Input validation
  if (!all(c("UB", "CB") %in% names(df))) {
    stop("Data frame must contain 'UB' and 'CB' columns")
  }

  df %>%
    distinct(CB, UB) %>%  # Get unique UMI-cell combinations
    group_by(CB) %>%
    summarise(UMI_per_CB = n(), .groups = "drop")
}

#' Filter mutations based on read count per UMI (internal function).
#'
#' This internal function filters mutations within a single UMI based on frequency threshold.
#' It keeps mutations that appear in at least freq_threshold proportion of reads for that UMI.
#' If no mutations pass the threshold and include_max is TRUE, it keeps the most abundant mutation.
#'
#' @param data A data frame of mutations for a single UMI.
#' @param freq_threshold Numeric threshold for mutation frequency (0-1).
#' @param include_max Logical indicating whether to include the maximum count mutation if none pass threshold.
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A filtered data frame of mutations.
#'
filter_mutations <- function(data, include_max = TRUE, freq_threshold) {
  # Filter mutations that meet the frequency threshold
  df <- data %>%
    filter(count > freq_threshold * read_count_per_UMI)

  # If no mutations pass threshold and include_max is TRUE, keep the most abundant
  if (include_max && nrow(df) == 0 && nrow(data) > 0) {
    max_count <- max(data$count)
    df <- data %>%
      filter(count == max_count)
  }

  return(df)
}

#' Filter mutations based on read count per UMI.
#'
#' This function filters mutations at the UMI level for single-cell data. For each UMI within each cell,
#' it keeps mutations that appear in at least freq_threshold proportion of reads supporting that UMI.
#' This helps remove PCR errors and sequencing noise.
#'
#' @param data A data frame of mutations with columns: CB, UB, type, start, length, mutate_to, count, read_count_per_UMI.
#' @param freq_threshold Numeric threshold for mutation frequency per UMI (0-1).
#' @param include_max Logical indicating whether to include the most abundant mutation if none pass threshold.
#'
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A filtered data frame of mutations.
#' @export
#'
filter_mutations_per_UMI <- function(data, include_max = TRUE, freq_threshold = 0.5) {
  # Input validation
  required_cols <- c("CB", "UB", "type", "start", "length", "mutate_to", "count", "read_count_per_UMI")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(freq_threshold) || freq_threshold < 0 || freq_threshold > 1) {
    stop("freq_threshold must be a number between 0 and 1")
  }

  # Aggregate counts per mutation per UMI
  df <- data %>%
    group_by(CB, UB, type, start, length, mutate_to) %>%
    summarise(
      count = sum(count),
      read_count_per_UMI = mean(read_count_per_UMI),
      .groups = "drop"
    )

  # Apply filtering per UMI
  df %>%
    group_by(CB, UB) %>%
    nest() %>%
    mutate(mutations = map(data, filter_mutations, include_max, freq_threshold)) %>%
    unnest(mutations) %>%
    select(-data) %>%
    ungroup()
}


#' Filter mutations based on UMI count per cell.
#'
#' This function filters mutations at the cell level for single-cell data. It keeps mutations
#' that are supported by at least freq_threshold proportion of UMIs within each cell.
#' This provides an additional layer of filtering to ensure mutations are not cell-level artifacts.
#'
#' @param data A data frame of mutations with columns: CB, type, start, length, mutate_to, UMI_per_CB.
#' @param freq_threshold Numeric threshold for mutation frequency per cell (0-1).
#'
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A filtered data frame of mutations.
#' @export
#'
filter_mutations_per_cell <- function(data, freq_threshold = 0.5) {
  # Input validation
  required_cols <- c("CB", "type", "start", "length", "mutate_to", "UMI_per_CB")
  if (!all(required_cols %in% names(data))) {
    stop("Data must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(freq_threshold) || freq_threshold < 0 || freq_threshold > 1) {
    stop("freq_threshold must be a number between 0 and 1")
  }

  data %>%
    group_by(CB, type, start, length, mutate_to) %>%
    summarise(
      UMI_count = n(),  # Number of UMIs supporting this mutation in this cell
      UMI_per_CB = mean(UMI_per_CB),
      .groups = "drop"
    ) %>%
    filter(UMI_count >= freq_threshold * UMI_per_CB)
}
