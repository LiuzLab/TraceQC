library(tidyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)


#' Display a circos plot with links for a given data frame.
#'
#' @param df a data frame that contains data to be visualized on the plot
#' @param title The main title of the plot
#' @param traceQC_input A TraceQC object
#' @param count_cutoff A cutoff to remove link whose log10-count are less than the value.
#'
#' @import circlize
#' @importFrom magrittr %>%
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Legend draw
#'
#' @return It doesn't generate any specific output.
#'
circular_chordgram <-
  function(df, title, ref, use_log_count=TRUE, count_cutoff = 1) {
    regions <- ref$regions
    target_start <- regions %>% filter(.data$region == "target") %>% pull(.data$start)
    target_end <- regions %>% filter(.data$region == "target") %>% pull(.data$end)
    refseq <-
      substr(ref$refseq, start = target_start, stop = target_end)

    col_fun <- colorRamp2(c(0,
                            ceiling(max(df$count))), c("yellow", "red"))
    df$color <- col_fun(df$count)
    l <- nchar(refseq)+1

    circos.par(start.degree = 90)
    circos.initialize(factors = factor(1), xlim = c(0, ceiling(l * 1.05)))
    circos.track(track.height = 0.1, ylim = c(0, 1))
    circos.axis(major.at = seq(0, l, 10), labels = seq(0, l, 10))
    df <-
      df %>% arrange(by = .data$count) %>% filter(.data$count >= count_cutoff)
    max_cnt <- ceiling(max(df$count))
    
    if (use_log_count) {
      for (i in 1:nrow(df)) {
        circos.link(
          1,
          df$start[i],
          1,
          df$end[i],
          h.ratio = 0.9,
          lwd = df$count[i] / 2,
          col = alpha(df$color[i], df$count[i] / max_cnt)
        )
      }} else {
      for (i in 1:nrow(df)) {
        circos.link(
          1,
          df$start[i],
          1,
          df$end[i],
          h.ratio = 0.9,
          lwd = 1,
          col = df$color[i]
        )
      }}

    region_names <- unique(regions$region)
    colors <- brewer.pal(length(region_names) + 1, "Set2")
    names(colors) <- region_names
    col = rep("black", nchar(refseq))
    for (i in 1:nrow(regions)) {
      col[(regions[i, "start"]):regions[i, "end"]] <- colors[regions[i,"region"]]
    }
    col <- col[target_start:target_end]

    circos.text(1:l,
                0.5,
                strsplit(refseq, split = "") %>% unlist(),
                col = col,
                cex = 1)


    title(title)
    circos.clear()
    lgd <-
      Legend(
        at = seq(floor(min(df$count)), ceiling(max(df$count))+1,
                 length.out = 5),
        col_fun = col_fun,
        title = "count"
      )
    draw(lgd, x = unit(0.15, "npc"), y = unit(0.15, "npc"))
  }

#' Display a circos plot with a histgoram for a given data frame.
#'
#' @param ref A traceQC object
#' @param title The main title of the plot.
#' @param df a data frame that contains data to be visualized on the plot.
#'
#' @importFrom magrittr %>%
#' @import circlize
#' @import dplyr
#' @importFrom ComplexHeatmap Legend draw
#' @importFrom grid gpar
#'
#' @return It doesn't generate any specific output.
#'
circular_histogram <- function(df, title, ref) {
  regions <- ref$regions
  target_start <- regions %>% filter(.data$region == "target") %>% pull(.data$start)
  target_end <- regions %>% filter(.data$region == "target") %>% pull(.data$end)
  refseq <-
    substr(ref$refseq, start = target_start, stop = target_end)

  scale <- df %>%
    group_by(.data$start) %>%
    summarise(count = sum(.data$count)) %>%
    ungroup()

  l <- nchar(refseq)+1

  circos.par(start.degree = 90)
  circos.initialize(factors = factor(1), xlim = c(0, ceiling(l * 1.05)))
  circos.track(
    track.height = 0.2,
    ylim = c(0, 1.2 * max(scale$count)),
    bg.border = NA
  )
  circos.track(track.height = 0.1, ylim = c(0, 1))
  circos.yaxis(track.index = 1, labels.cex	= 0.5)
  circos.axis(
    h = "bottom",
    labels.facing = "reverse.clockwise",
    direction = "inside",
    major.at = seq(0, l, 10),
    labels = seq(0, l, 10)
  )
  region_names <- unique(regions$region)
  colors <- brewer.pal(length(region_names) + 1, "Set2")
  names(colors) <- region_names
  col = rep("black", nchar(refseq))
  for (i in 1:nrow(regions)) {
    col[(regions[i, "start"]):regions[i, "end"]] <- colors[regions[i,"region"]]
  }
  col <- col[target_start:target_end]
  circos.text(1:l,
              0.5,
              strsplit(refseq, split = "") %>% unlist(),
              col = col,
              cex = 1)

  colors <- c("#d7191c", "#fdae61", "#abd9e9", "#2c7bb6")
  names(colors) <- c("A", "C", "G", "T")

  for (i in 1:nrow(df)) {
    circos.rect(
      df$start[i] - 0.4,
      df$y[i] - df$count[i] + 1,
      df$start[i] + 0.4,
      df$y[i],
      col = colors[df$mutate_to[i]],
      border = NA,
      track.index = 1
    )
  }

  title(title)
  circos.clear()

  lgd <- Legend(
    at = c("A", "C", "G", "T"),
    type = "points",
    legend_gp = gpar(col = c(
      "#d7191c", "#fdae61", "#abd9e9", "#2c7bb6"
    )),
    title_position = "topleft",
    title = "Nucleotide"
  )

  draw(lgd, x = unit(0.15, "npc"), y = unit(0.15, "npc"))

}

#' Display a circos plot that shows overall deletion pattern across the barcodes.
#'
#' @param mutations A mutations dataframe.
#' @param ref A reference object.
#' @param count_cutoff A cutoff to remove link whose log10-count are less than the value.
#'
#' @importFrom magrittr %>%
#' @import circlize
#' @import dplyr
#'
#' @return It doesn't generate any specific output.
#'
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_deletion_hotspot(example_obj)
#'
plot_deletion_hotspot <- function(mutations, ref, use_log_count = TRUE, count_cutoff = 1) {
  if (!("count" %in% names(mutations))) {mutations$count <- 1}
  deletions <- mutations %>%
    filter(.data$type == "deletion") %>%
    group_by(.data$start, .data$length) %>%
    summarise(count = sum(.data$count)) %>%
    mutate(end = .data$start + .data$length) %>%
    ungroup %>%
    mutate(id = 1:n()) %>%
    mutate(log10_count = log10(.data$count))
  
  if (use_log_count) {
    deletions$count = log10(deletions$count)
  }
  
  circular_chordgram(df = deletions,
                     title = "Deletions",
                     ref,
                     use_log_count = use_log_count,
                     count_cutoff)
}

#' Display a circos plot that shows overall insertion pattern across the barcodes.
#'
#' @param mutations A mutations dataframe.
#' @param ref A reference object.
#' @param count_cutoff A cutoff to remove link whose log10-count are less than the value.
#'
#' @importFrom magrittr %>%
#' @import circlize
#' @import dplyr
#'
#' @return It won't return any specific object.
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_insertion_hotspot(example_obj)
#'
plot_insertion_hotspot <-
  function(mutations, ref, use_log_count=TRUE, count_cutoff = 1) {
    if (!("count" %in% names(mutations))) {mutations$count <- 1}
    insertions <- mutations %>%
      filter(.data$type == "insertion") %>%
      group_by(.data$start, .data$length) %>%
      summarise(count = sum(.data$count)) %>%
      mutate(end = .data$start + .data$length) %>%
      ungroup %>%
      mutate(id = 1:n())
    if (use_log_count) {
      insertions$count = log10(insertions$count)
    }

    circular_chordgram(df = insertions,
                       title = "Insertions",
                       ref,
                       use_log_count = use_log_count,
                       count_cutoff)
  }

#' Display a mutation hotspot circos plot.
#'
#' The circos plot shows the frequency of mutation events for each nucleotide.
#'
#' @param mutations A mutations dataframe
#' @param ref A reference object
#'
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom purrr map_chr
#'
#' @return It won't return any specific object.
#'
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_point_mutation_hotspot(example_obj)
#'
plot_point_substitution_hotspot <- function(mutations,ref) {
  if (!("count" %in% names(mutations))) {mutations$count <- 1}
  substitutions <- mutations %>% filter(.data$type == "substitution") %>%
    group_by(.data$start, .data$length, .data$mutate_to) %>%
    summarise(count = sum(.data$count)) %>%
    ungroup %>%
    group_by(.data$start) %>%
    mutate(mutate_to = as.character(.data$mutate_to)) %>%
    arrange(.data$count) %>%
    mutate(y = cumsum(.data$count)) %>%
    ungroup
  circular_histogram(substitutions, "Substitutions", ref)
}

