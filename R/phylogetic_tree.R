#' Convert A TraceQC object to a matrix
#'
#' @param TraceQC_input A TraceQC object
#' @param count_cutoff The barcode cutoff. The function only considers barcodes
#' whose raw counts are equal or greater than the cutoff.
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @return It returns a list that contains a character matrix, a sequence
#' information of each barcode, and a data frame that contains code of each
#' barcode sequence.
#'
#' @export
#'
build_binary_table <- function(TraceQC_input,
                               count_cutoff = 5) {
  df <- TraceQC_input$mutation %>%
    dplyr::filter(.data$count>count_cutoff,.data$type!="unmutated")
  seq_id <- df %>% group_by(.data$target_seq, .data$count) %>%
    summarise() %>%
    ungroup %>%
    mutate(sequenceID = 1:n())

  unique_events <- df %>%
    group_by(.data$type, .data$start, .data$length, .data$mutate_to) %>%
    summarise() %>%
    ungroup %>%
    mutate(mutationID = 1:n())

  event_str <- unique_events %>%
    mutate(type = ifelse(.data$type == "deletion", "D",
                         ifelse(.data$type == "insertion", "I", "M"))) %>%
    mutate(mutate_to = ifelse(.data$type == "D", "", .data$mutate_to) ) %>%
    mutate(length = ifelse(.data$type == "D", .data$length, "")) %>%
    mutate(event = ifelse(.data$type == "D", .data$length, .data$mutate_to)) %>%
    mutate(code = paste(.data$type, .data$start, "[", .data$event, "]", sep = "")) %>%
    select(.data$mutationID, .data$code)

  nr <- nrow(seq_id)
  nc <- nrow(unique_events)

  df <- left_join(df, seq_id) %>%
    left_join(unique_events) %>%
    mutate(coord = nr * (.data$mutationID - 1) + .data$sequenceID)

  binary_table <- matrix(data = 0,
                            nrow = nr,
                            ncol = nc)
  binary_table[df$coord] <- 1
  rownames(binary_table) <- 1:nr


  return(list(
    binary_table = binary_table,
    seq_id = seq_id,
    event_str = event_str
    ))
}


#' Convert a binary table of mutation events to a phyDat object
#'
#' @param tree_input A variable which was created by `build_binary_table'
#' @import dplyr
#' @import phangorn
#' @importFrom magrittr %>%
#' @import treeio
#' @import ggtree
#'
#' @return
#' @export
#'
table_to_phyDat <- function(tree_input) {
  mat <- tree_input$binary_table
  colnames(mat) <- tree_input$event_str$code

  get_events <- function(mat, i) {
    paste(names(which(mat[i,]==1)), collapse = ",")
  }

  rownames(mat) <-
    sapply(1:nrow(mat), function(i) get_events(mat, i))

  get_events_only_types <- function(x) {
    if(length(x) == 1) {
      if(is.na(x)) {
        return(NA)
      }
      ret <- c()
      if(stringr::str_detect(x, "M")) ret <- c(ret, "Mutation")
      if(stringr::str_detect(x, "I")) ret <- c(ret, "Insertion")
      if(stringr::str_detect(x, "D")) ret <- c(ret, "Deletion")
      paste(ret, collapse = "+")
    } else {
      sapply(x, function(i) get_events_only_types(i))
    }
  }

  data <- phyDat(data=mat,type="USER",levels=c(0,1))
  dm <- dist.hamming(data)
  tree_UPGMA <- upgma(dm)
  # tree_pars <- optim.parsimony(tree_UPGMA, data)

  # df_cnt <- tibble::tibble(code = rownames(mat),
  #                          cnt = tree_input$seq_id$count)

#   tree_obj <- tibble::as_tibble(tree_UPGMA) %>%
#     filter(!is.na(label)) %>%
#     dplyr::mutate(type = get_events_only_types(.data$label)) %>%
#     dplyr::left_join(df_cnt, by=c("label"="code")) %>%
#     treeio::as.treedata()
#
#   tree_obj
  tree_UPGMA
}

#' Plot a phylogenetic tree of a tree object using ggtree
#'
#' @param tree_obj A variable which was created by `table_to_phyDat'
#'
#' @return It returns a ggplot2 object of a phylogenetic tree.
#' @export
#'
#' @examples
#' library(TraceQC)
#' library(magrittr)
#' library(cowplot)
#' library(ggplot2)
#'
#' data("example_obj_list")
#' p1 <- example_obj_list$day00 %>%
#' build_binary_table %>%
#' table_to_phyDat %>%
#' plot_phylogenetic_tree
#' p2 <- example_obj_list$day02 %>%
#' build_binary_table %>%
#' table_to_phyDat %>%
#' plot_phylogenetic_tree
#' p3 <- example_obj_list$day14 %>%
#' build_binary_table %>%
#' table_to_phyDat %>%
#' plot_phylogenetic_tree
#'
#' plot_grid(p1 + ggtitle("Day 0"),
#'           p2 + ggtitle("Day 2"),
#'           p3 + ggtitle("Day 14"),
#'           nrow=1)
#'
plot_phylogenetic_tree <- function(tree_obj) {
  span <- max(ape::node.depth.edgelength(tree_obj@phylo))
  p <- tree_obj %>%
    ggtree::ggtree(layout="rectangular") +
    ggplot2::xlim(0, span * 1.2) +
    ggtree::geom_tiplab(size=2,
                        ggplot2::aes_string(fill="type"),
                        #offset=0.1,
                        geom = "label") +
    ggtree::geom_tippoint(ggplot2::aes_string(color = "log10(cnt)"), size=1) +
    ggplot2::scale_color_gradient(low = "yellow", high = "red", na.value = NA) +
    ggtree::theme_tree2()
  p
}
