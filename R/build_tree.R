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
  df <- filter(TraceQC_input$mutation,count>5,type!="unmutated")
  seq_id <- group_by(df, target_seq, count) %>%
    summarise() %>%
    ungroup %>%
    mutate(sequenceID = 1:n())

  unique_events <- df %>%
    group_by(type, start, length, mutate_to) %>%
    summarise() %>%
    ungroup %>%
    mutate(mutationID = 1:n())

  nr <- nrow(seq_id)
  nc <- nrow(unique_events)

  df <- left_join(df, seq_id) %>%
    left_join(unique_events) %>%
    mutate(coord = nr * (mutationID - 1) + sequenceID)

  character_table <- matrix(data = 0,
                            nrow = nr,
                            ncol = nc)
  character_table[df$coord] <- 1
  rownames(character_table) <- 1:nr
  return(list(character_table = character_table,
              seq_id = seq_id))
}
