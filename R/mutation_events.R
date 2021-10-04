#' Creating a data frame of mutation events.
#'
#' @param insertions A list that contains insertion events.
#' @param deletions A list that contains deletion events.
#' @param mutations A list that contains mutation (substitution) events.
#' @param target_seq A list that contains the alignment for each event.
#' @param score A list that contains the alignment score for each event.
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
           score,
           read_count) {
    insertions <- insertions %>%
      as.data.frame() %>%
      mutate(target_seq = target_seq) %>%
      mutate(
        length = .data$end - .data$start + 1,
        type = "insertion",
        mutate_to = substr(target_seq, start = .data$start, stop = .data$end)
      )
    deletions <- deletions %>%
      as.data.frame() %>%
      mutate(target_seq = target_seq) %>%
      mutate(length = .data$end - .data$start + 1,
             type = "deletion",
             mutate_to = "-")
    mutations <- data.frame(start = mutations, end = mutations) %>%
      mutate(target_seq = target_seq) %>%
      mutate(
        length = 1,
        type = "substitution",
        mutate_to = substr(target_seq, start = .data$start, stop = .data$end)
      )
    return (
      rbind(insertions, deletions, mutations) %>%
        arrange(.data$start) %>%
        mutate(tmp = lag(as.integer(.data$type == "insertion") * length)) %>%
        mutate(tmp = replace_na(.data$tmp, 0)) %>%
        mutate(align = cumsum(.data$tmp)) %>%
        mutate(
          start = .data$start - .data$align,
          target_seq = .data$target_seq,
          alignment_score = score,
          count = read_count
        ) %>%
        select(target_seq, alignment_score, type, start, length, mutate_to, count)
    )
  }

#' Identifying mutation events.
#'
#' @param traceQC_input A TraceQC object.
#' @param use_CPM Use count per million
#' @param alignment_score_cutoff Minimum cutoff for alignment score
#' @param abundance_cutoff Minimum cutoff for read count. This parameter are used with use_CPM.
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
seq_to_character <- function(aligned_reads,
                             use_CPM,
                             alignment_score_cutoff = 0,
                             abundance_cutoff = 0) {
  if (!"count" %in% names(aligned_reads)) {
    aligned_reads$count <- 1
    abundance_cutoff <- 0}
  
  if(use_CPM) {
    aligned_reads$count <- aligned_reads$count * 1e6 / sum(aligned_reads$count)
    abundance_cutoff <- 1e6 * abundance_cutoff
  } else {
    abundance_cutoff <- sum(aligned_reads$count) * abundance_cutoff
  }
  unmutated <- aligned_reads %>%
    filter(.data$target_seq == .data$target_ref)
  aligned_reads <- filter(aligned_reads,count>abundance_cutoff,
                          score>alignment_score_cutoff)

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
  scores <- aligned_reads$score

  mutation_df <- NULL

  mutation_df <- mapply(
    find_position,
    all_insertions,
    all_deletions,
    all_mutations,
    target_seqs,
    scores,
    read_counts,
    SIMPLIFY = FALSE
  ) %>% bind_rows()

  if (nrow(unmutated)>0) {
  mutation_df <- rbind(
    data.frame(
      target_seq = unmutated$target_seq,
      alignment_score = unmutated$score,
      type = "unmutated",
      start = 0,
      length = 0,
      mutate_to = "-",
      count = unmutated$count
    ),
    mutation_df
  )}

  return(mutation_df)
}


#' format a mutation data frame for output.
#' 
#' @param mutations A data frame of mutations. The output of seq_to_character.
#' 
#' @importFrom tidyr nest
#' @importFrom purrr map_chr
#' @importFrom purrr map_dbl
#' @return A formatted data frame of mutations.
#' @export
#'
format_mutation_df <- function(mutation_df, is_singlecell) {
  if (is_singlecell) {
    out_df <- mutate(mutation_df,alt=case_when(mutate_to=="-" ~ "",
                                   TRUE ~ mutate_to)) %>%
      # filter(type!="unmutated") %>%
      group_by(type,start,length,alt) %>%
      nest() %>%
      mutate(cell=map_chr(data,function(x) {paste(unique(x$CB),collapse=",")})) %>%
      ungroup %>%
      select(-data) %>%
      mutate(character=case_when(type=="deletion"|type=="unmutated" ~ sprintf("%s:S%sL%s",toupper(
                        substr(type,start=1,stop=1)),start,length),
                                 TRUE ~ sprintf("%s:S%sL%s->%s",toupper(
                                   substr(type,start=1,stop=1)),start,length,alt)))
  } else {
    out_df <- mutate(mutation_df,alt=case_when(mutate_to=="-" ~ "",
                                               TRUE ~ mutate_to)) %>%
      group_by(type,start,length,alt) %>%
      summarise(count=sum(count)) %>%
      ungroup %>%
      mutate(character=case_when(type=="deletion"|type=="unmutated" ~ sprintf("%s:S%sL%s",toupper(
        substr(type,start=1,stop=1)),start,length),
        TRUE ~ sprintf("%s:S%sL%s->%s",toupper(
          substr(type,start=1,stop=1)),start,length,alt)))}
  return(relocate(out_df,character,type,start,length,alt))
}

