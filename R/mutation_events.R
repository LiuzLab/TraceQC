#' Creating a data frame of mutation events for a single aligned read.
#'
#' This internal function processes the raw mutation data (insertions, deletions, substitutions)
#' for a single aligned read and converts it into a standardized data frame format.
#' It adjusts mutation positions to account for insertions that shift downstream positions.
#'
#' @param insertions A matrix of start/end positions for insertion events (from str_locate_all).
#' @param deletions A matrix of start/end positions for deletion events (from str_locate_all).
#' @param mutations A vector of positions for substitution events.
#' @param target_seq The aligned target sequence string.
#' @param score The alignment score for this read.
#' @param read_count The count/abundance of this read.
#'
#' @import dplyr
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#'
#' @return A data frame with columns: target_seq, alignment_score, type, start, length, mutate_to, count.
#'
find_position <- function(insertions, deletions, mutations, target_seq, score, read_count) {
  # Process insertions: extract sequences and create data frame
  insertions_df <- if (nrow(insertions) > 0) {
    insertions %>%
      as.data.frame() %>%
      mutate(
        target_seq = target_seq,
        length = .data$end - .data$start + 1,
        type = "insertion",
        mutate_to = substr(target_seq, start = .data$start, stop = .data$end)
      )
  } else {
    data.frame(start = integer(), end = integer(), target_seq = character(),
               length = integer(), type = character(), mutate_to = character())
  }

  # Process deletions: mark as gaps
  deletions_df <- if (nrow(deletions) > 0) {
    deletions %>%
      as.data.frame() %>%
      mutate(
        target_seq = target_seq,
        length = .data$end - .data$start + 1,
        type = "deletion",
        mutate_to = "-"
      )
  } else {
    data.frame(start = integer(), end = integer(), target_seq = character(),
               length = integer(), type = character(), mutate_to = character())
  }

  # Process substitutions: single base changes
  mutations_df <- if (length(mutations) > 0) {
    data.frame(start = mutations, end = mutations) %>%
      mutate(
        target_seq = target_seq,
        length = 1L,
        type = "substitution",
        mutate_to = substr(target_seq, start = .data$start, stop = .data$end)
      )
  } else {
    data.frame(start = integer(), end = integer(), target_seq = character(),
               length = integer(), type = character(), mutate_to = character())
  }

  # Combine all mutation types
  combined <- bind_rows(insertions_df, deletions_df, mutations_df)

  if (nrow(combined) == 0) {
    return(data.frame(target_seq = character(), alignment_score = numeric(),
                      type = character(), start = integer(), length = integer(),
                      mutate_to = character(), count = numeric()))
  }

  # Adjust positions: insertions shift downstream mutations
  combined <- combined %>%
    arrange(.data$start) %>%
    mutate(
      # Calculate cumulative shift from preceding insertions
      tmp = lag(as.integer(.data$type == "insertion") * .data$length),
      tmp = replace_na(.data$tmp, 0),
      align = cumsum(.data$tmp),
      # Adjust start position
      start = .data$start - .data$align,
      alignment_score = score,
      count = read_count
    ) %>%
    select(target_seq, alignment_score, type, start, length, mutate_to, count)

  return(combined)
}

#' Identifying mutation events from aligned reads.
#'
#' This function analyzes aligned sequencing reads to identify and characterize mutations
#' (insertions, deletions, substitutions) relative to the reference sequence. It filters
#' reads based on alignment score and abundance, then processes each read to extract
#' mutation information.
#'
#' @param aligned_reads A data frame of aligned reads with columns: target_seq, target_ref, score, count (optional).
#' @param use_CPM Logical indicating whether to use counts per million normalization.
#' @param alignment_score_cutoff Numeric minimum alignment score to keep a read.
#' @param abundance_cutoff Numeric minimum abundance cutoff (absolute or CPM depending on use_CPM).
#'
#' @importFrom parallel mcmapply
#' @importFrom stringr str_locate_all
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item `target_seq`: The aligned sequence.
#'     \item `alignment_score`: The alignment score.
#'     \item `type`: Mutation type ("insertion", "deletion", "substitution", "unmutated").
#'     \item `start`: Starting position of the mutation (1-based).
#'     \item `length`: Length of the mutation.
#'     \item `mutate_to`: The mutated sequence (or "-" for deletions).
#'     \item `count`: Read count/abundance.
#'   }
#' @export
#' @examples
#' # Assuming aligned_reads is loaded from sequence_alignment output
#' unique_sequences <- group_by(aligned_reads, target_seq, target_ref) %>%
#'   summarise(count = n(), score = max(score)) %>%
#'   ungroup()
#' alignment_score_cutoff <- 50  # Example threshold
#' mutations <- seq_to_character(unique_sequences, use_CPM = FALSE,
#'                               alignment_score_cutoff = alignment_score_cutoff)
#'
seq_to_character <- function(aligned_reads,
                             use_CPM = FALSE,
                             alignment_score_cutoff = 0,
                             abundance_cutoff = 0) {
  # Input validation
  required_cols <- c("target_seq", "target_ref", "score")
  if (!all(required_cols %in% names(aligned_reads))) {
    stop("aligned_reads must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!is.numeric(alignment_score_cutoff) || alignment_score_cutoff < 0) {
    stop("alignment_score_cutoff must be a non-negative number")
  }
  if (!is.numeric(abundance_cutoff) || abundance_cutoff < 0) {
    stop("abundance_cutoff must be a non-negative number")
  }

  # Ensure count column exists
  if (!"count" %in% names(aligned_reads)) {
    aligned_reads$count <- 1
    abundance_cutoff <- 0
  }

  # Apply CPM normalization if requested
  if (use_CPM) {
    total_reads <- sum(aligned_reads$count)
    if (total_reads == 0) stop("Total read count is zero, cannot compute CPM")
    aligned_reads$count <- aligned_reads$count * 1e6 / total_reads
    abundance_cutoff <- 1e6 * abundance_cutoff
  }

  # Separate unmutated reads (perfect matches)
  unmutated <- aligned_reads %>%
    filter(.data$target_seq == .data$target_ref)

  # Filter reads based on abundance and alignment score
  aligned_reads <- aligned_reads %>%
    filter(.data$count > abundance_cutoff, .data$score > alignment_score_cutoff)

  # Early return if no reads pass filters
  if (nrow(aligned_reads) == 0) {
    warning("No reads passed the filtering criteria")
    return(data.frame(
      target_seq = character(),
      alignment_score = numeric(),
      type = character(),
      start = integer(),
      length = integer(),
      mutate_to = character(),
      count = numeric()
    ))
  }

  # Identify mutation positions for each read
  # Insertions: gaps in reference (- in target_ref)
  all_insertions <- str_locate_all(aligned_reads$target_ref, "-+")
  # Deletions: gaps in target sequence (- in target_seq)
  all_deletions <- str_locate_all(aligned_reads$target_seq, "-+")
  # Substitutions: mismatched non-gap positions
  all_mutations <- mapply(
    function(ref, seq) {
      which(ref != seq & ref != "-" & seq != "-" & seq != "N")
    },
    strsplit(aligned_reads$target_ref, ""),
    strsplit(aligned_reads$target_seq, ""),
    SIMPLIFY = FALSE
  )

  # Extract read metadata
  read_counts <- aligned_reads$count
  target_seqs <- aligned_reads$target_seq
  scores <- aligned_reads$score

  # Process mutations for each read in parallel
  mutation_df <- mcmapply(
    find_position,
    all_insertions,
    all_deletions,
    all_mutations,
    target_seqs,
    scores,
    read_counts,
    SIMPLIFY = FALSE
  ) %>% bind_rows()

  # Add unmutated reads if any
  if (nrow(unmutated) > 0) {
    unmutated_df <- data.frame(
      target_seq = unmutated$target_seq,
      alignment_score = unmutated$score,
      type = "unmutated",
      start = 0L,
      length = 0L,
      mutate_to = "-",
      count = unmutated$count
    )
    mutation_df <- bind_rows(unmutated_df, mutation_df)
  }

  return(mutation_df)
}


#' Format a mutation data frame for output.
#'
#' This function aggregates and formats mutation data for final output, creating
#' standardized mutation identifiers and optionally including cell barcode information
#' for single-cell data.
#'
#' @param mutation_df A data frame of mutations, typically the output of seq_to_character.
#' @param is_singlecell Logical indicating if the data is from single-cell sequencing.
#'   If TRUE, expects a 'CB' column for cell barcodes.
#'
#' @importFrom tidyr nest
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @importFrom magrittr %>%
#' @import dplyr
#' @return A formatted data frame with columns: character, type, start, length, alt, [cell], [count].
#'   The 'character' column contains standardized mutation identifiers like "I:S10L3->ATG".
#' @export
#' @examples
#' # For bulk data
#' out_df <- format_mutation_df(mutations, is_singlecell = FALSE)
#' # For single-cell data (assuming mutations has CB column)
#' out_df <- format_mutation_df(mutations, is_singlecell = TRUE)
#'
format_mutation_df <- function(mutation_df, is_singlecell = FALSE) {
  # Input validation
  required_cols <- c("type", "start", "length", "mutate_to")
  if (!all(required_cols %in% names(mutation_df))) {
    stop("mutation_df must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (is_singlecell && !"CB" %in% names(mutation_df)) {
    stop("For single-cell data, mutation_df must contain a 'CB' column")
  }

  # Create standardized mutation identifiers
  out_df <- mutation_df %>%
    mutate(
      alt = case_when(
        mutate_to == "-" ~ "",
        TRUE ~ mutate_to
      )
    )

  if (is_singlecell) {
    # Group by mutation characteristics and aggregate cell information
    out_df <- out_df %>%
      group_by(type, start, length, alt) %>%
      nest() %>%
      mutate(
        cell = map_chr(data, ~ paste(unique(.x$CB), collapse = ","))
      ) %>%
      ungroup() %>%
      select(-data)
  } else {
    # For bulk data, sum counts across identical mutations
    out_df <- out_df %>%
      group_by(type, start, length, alt) %>%
      summarise(count = sum(count), .groups = "drop")
  }

  # Create character representation of mutations
  out_df <- out_df %>%
    mutate(
      character = case_when(
        type == "deletion" | type == "unmutated" ~
          sprintf("%s:S%sL%s", toupper(substr(type, 1, 1)), start, length),
        TRUE ~
          sprintf("%s:S%sL%s->%s", toupper(substr(type, 1, 1)), start, length, alt)
      )
    ) %>%
    relocate(character, type, start, length, alt)

  return(out_df)
}

