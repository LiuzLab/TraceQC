#' Convert A TraceQC object to a matrix
#'
#' @param TraceQC_input A TraceQC object
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @return It returns a list that contains a character matrix and a sequence
#' information.
#'
#' @export
#'
build_character_table <- function(TraceQC_input) {
  df <- TraceQC_input$mutation %>% filter(.data$count>5,.data$type!="unmutated")
  seq_id <- df %>% group_by(.data$target_seq, .data$count) %>%
    summarise() %>%
    ungroup %>%
    mutate(sequenceID = 1:n())

  unique_events <- df %>%
    group_by(.data$type, .data$start, .data$length, .data$mutate_to) %>%
    summarise() %>%
    ungroup %>%
    mutate(mutationID = 1:n())

  nr <- nrow(seq_id)
  nc <- nrow(unique_events)

  df <- left_join(df, seq_id) %>%
    left_join(unique_events) %>%
    mutate(coord = nr * (.data$mutationID - 1) + .data$sequenceID)

  character_table <- matrix(data = 0,
                            nrow = nr,
                            ncol = nc)
  character_table[df$coord] <- 1
  rownames(character_table) <- 1:nr
  return(list(character_table = character_table,
              seq_id = seq_id))
}
