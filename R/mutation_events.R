#' Creating a data frame of mutation events.
#'
#' @param insertions A list that contains insertion events.
#' @param deletions A list that contains deletion events.
#' @param mutations A list that contains mutation (substitution) events.
#' @param target_seq A list that contains the alignment for each event.
#' @param read_count A vector that contains counts for each event.
#'
#' @import dplyr
#' @importFrom tidyr replace_na
#' @importFrom magrittr %>%
#'
#' @return A data frame that contains the event information.
#'
find_position <-
  function(insertions,
           deletions,
           mutations,
           target_seq,
           read_count) {
    insertions <- insertions %>%
      as.data.frame() %>%
      mutate(target_seq = target_seq) %>%
      mutate(
        length = end - start + 1,
        type = "insertion",
        mutate_to = substr(target_seq, start = start, stop = end)
      )
    deletions <- deletions %>%
      as.data.frame() %>%
      mutate(target_seq = target_seq) %>%
      mutate(length = end - start + 1,
             type = "deletion",
             mutate_to = "-")
    mutations <- data.frame(start = mutations, end = mutations) %>%
      mutate(target_seq = target_seq) %>%
      mutate(
        length = 1,
        type = "mutation",
        mutate_to = substr(target_seq, start = start, stop = end)
      )
    return (
      rbind(insertions, deletions, mutations) %>%
        arrange(start) %>%
        mutate(tmp = lag(as.integer(
          type == "insertion"
        ) * length)) %>%
        mutate(tmp = replace_na(tmp, 0)) %>%
        mutate(align = cumsum(tmp)) %>%
        mutate(
          start = start - align,
          target_seq = target_seq,
          count = read_count
        ) %>%
        select(target_seq, type, start, length, mutate_to, count)
    )
  }

#' Identifying mutation events.
#'
#' @param traceQC_input A TraceQC object.
#' @param ncores The number of cores for parallelization.
#'
#' @importFrom parallel mcmapply
#' @importFrom stringr str_locate_all
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A data frame that contains the columns:
#' \itemize{
#'   \item `type': The type of mutation.
#'   \item `start': The starting position of mutation event.
#'   \item `length': The length of mutation.
#'   \item `mutation_to': A string that shows what mutation is occurred.
#'   \item `count': The total number of the mutation events from the sample.
#' }
#' @export
#'
seq_to_character <- function(traceQC_input,
                             ncores = 4) {
  aligned_reads <- traceQC_input$aligned_reads
  col_names <- names(aligned_reads)
  col_names <-
    col_names[!(col_names %in% c("name", "seq", "ref", "score"))]

  aligned_reads <- aligned_reads %>%
    group_by_at(col_names) %>%
    summarise(count = n()) %>%
    ungroup

  all_insertions <- str_locate_all(aligned_reads$target_ref, "-+")
  all_deletions <- str_locate_all(aligned_reads$target_seq, "-+")
  all_mutations <-
    mapply(
      function(x, y) {
        which(x != y & x != "-" & y != "-" & y != "N")
      },
      strsplit(aligned_reads$target_ref, ""),
      strsplit(aligned_reads$target_seq, "")
    )
  read_counts <- aligned_reads$count
  target_seqs <- aligned_reads$target_seq

  mutation_df <- mcmapply(
    find_position,
    all_insertions,
    all_deletions,
    all_mutations,
    target_seqs,
    read_counts,
    SIMPLIFY = FALSE,
    mc.cores = ncores
  ) %>% bind_rows()

  unmutated <- filter(aligned_reads, target_seq == target_ref)
  mutation_df <- rbind(
    data.frame(
      target_seq = unmutated$target_seq,
      type = "unmutated",
      start = 0,
      length = 0,
      mutate_to = "-",
      count = unmutated$count
    ),
    mutation_df
  )

  return(mutation_df)
}

