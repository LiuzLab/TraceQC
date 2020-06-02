
#' Create a list of TraceQC objects from a list
#'
#' @param samples A list that contains FASTQ file path.
#' @param ref The path to the reference file.
#' @param ncores The number of cores for parallel processing.
#'
#' @return It returns a list of TraceQC object
#' @export
#'
#' @examples
#' \dontrun{
#' samples <- list(
#'   "day00" = system.file("extdata", "test_data", "fastq",
#'                         "example_0d.fastq", package="TraceQC"),
#'
#'   "day02" = system.file("extdata", "test_data", "fastq",
#'                         "example_2d.fastq", package="TraceQC"),
#'
#'   "day14" = system.file("extdata", "test_data", "fastq",
#'                         "example_14d.fastq", package="TraceQC")
#' )
#'
#' ref <- system.file("extdata", "test_data", "ref",
#'                    "ref.txt", package="TraceQC")
#'
#' obj_list <- create_obj_list(samples, ref)
#' }
create_obj_list <- function(samples, ref, ncores = 4) {
  obj_list <- list()
  for(name in names(samples)) {
    path <- samples[[name]]
    obj <- TraceQC(path, ref, ncores= ncores)
    obj_list[[name]] <- obj
  }
  obj_list
}



#' Plot the mutation pattern of given samples using a line plot.
#'
#' @param obj_list A list created by `create_obj_list'.
#'
#' @return It returns a ggplot object that contains line plot.
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
#' data(example_obj_list)
#' plot_mutation_pattern_lineplot(example_obj_list)
plot_mutation_pattern_lineplot <- function(obj_list) {
  df <- tibble()
  for (l in names(obj_list)) {
    df <- bind_rows(
      df,
      obj_list[[l]]$mutation %>% group_by(type) %>%
        summarize(freq = sum(count)) %>% ungroup() %>%
        mutate(freq = freq / sum(freq) * 100) %>%
        mutate(label = l)
    )
  }
  df$label <- factor(df$label, levels=names(obj_list))

  ggplot(df, aes(x = label, y = freq)) +
    geom_point(aes(color = type, shape = type), size = 4) +
    geom_path(aes(group = type, color = type)) +
    ylab("Percentage (%)") +
    theme_classic()
}


#' Plot the mutation pattern of given samples using a violin plot.
#'
#' @param obj_list A list created by `create_obj_list'.
#'
#' @return It returns a ggplot object that contains violin plot.
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
#' data(example_obj_list)
#' plot_mutation_pattern_violinplot(example_obj_list)
plot_mutation_pattern_violinplot <- function(obj_list) {
  df <- tibble()
  for (l in names(obj_list)) {
    df <- bind_rows(
      df,
      obj_list[[l]]$mutation %>%
        filter(type != "unmutated") %>%
        mutate(freq = count / sum(count) * 100) %>%
        mutate(label = l)
    )
  }

  df$label <- factor(df$label, levels=names(obj_list))
  df %>%
    ggplot(aes(x = label, y = freq)) +
    geom_violin(aes(fill = label)) +
    ylab("Percentage (%)") +
    facet_grid(type ~ ., scales = "free_y") +
    theme_classic()
}
