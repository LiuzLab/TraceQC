#' Title
#'
#' @param df
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples
build_character_table <- function(df) {
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

#' Function to locate where the Quality Control file is stored.
#'
#' @param input_file A path of FASTQ file.
#' @param qc_dir The directory path where Quality Control files are.
#'
#' @return A path to the corresponded zipped Quality Control file.
#' @export
get_qcpath <- function(input_file, qc_dir) {
  qc_path <- sub(".fastq", "_fastqc.zip", basename(input_file))
  file.path(qc_dir, qc_path)
}
