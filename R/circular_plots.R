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
  function(df, title, traceQC_input, count_cutoff = 1) {
    regions <- traceQC_input$regions
    target_start <- filter(regions, region == "target") %>% pull(start)
    target_end <- filter(regions, region == "target") %>% pull(end)
    refseq <-
      substr(traceQC_input$refseq, start = target_start + 1, stop = target_end)

    col_fun <- colorRamp2(c(0,
                            ceiling(max(df$log10_count))), c("yellow", "red"))
    df$color <- col_fun(df$log10_count)
    l <- nchar(refseq) + 1

    circos.par(start.degree = 90)
    circos.initialize(factors = factor(1), xlim = c(0, ceiling(l * 1.05)))
    circos.track(track.height = 0.1, ylim = c(0, 1))
    circos.axis(major.at = seq(0, l, 10), labels = seq(0, l, 10))
    df <-
      df %>% arrange(by = log10_count) %>% filter(log10_count >= count_cutoff)
    max_cnt <- ceiling(max(df$log10_count))
    for (i in 1:nrow(df)) {
      circos.link(
        1,
        df$start[i],
        1,
        df$end[i],
        h.ratio = 0.9,
        lwd = df$log10_count[i] / 2,
        col = alpha(df$color[i], df$log10_count[i] / max_cnt)
      )
    }

    colors <- brewer.pal(nrow(regions) + 1, "Set2")
    col = rep("black", nchar(refseq))
    for (i in 1:nrow(regions)) {
      col[(regions[i, "start"] + 1):regions[i, "end"]] <- colors[i]
    }
    col <- col[target_start + 1:target_end]

    circos.text(1:l,
                0.5,
                strsplit(refseq, split = "") %>% unlist(),
                col = col,
                cex = 1)


    title(title)
    circos.clear()
    lgd <-
      Legend(
        at = seq(floor(min(df$log10_count)), ceiling(max(df$log10_count)),
                 length.out = 5),
        col_fun = col_fun,
        title = "log 10 count"
      )
    draw(lgd, x = unit(0.15, "npc"), y = unit(0.15, "npc"))
  }

#' Display a circos plot with a histgoram for a given data frame.
#'
#' @param traceQC_input A traceQC object
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
circular_histogram <- function(df, title, traceQC_input) {
  regions <- traceQC_input$regions
  target_start <- filter(regions, region == "target") %>% pull(start)
  target_end <- filter(regions, region == "target") %>% pull(end)
  refseq <-
    substr(traceQC_input$refseq, start = target_start + 1, stop = target_end)

  scale <- df %>%
    group_by(start) %>%
    summarise(count = sum(count)) %>%
    ungroup

  l <- nchar(refseq) + 1

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
  colors <- brewer.pal(nrow(regions) + 1, "Set2")
  col = rep("black", nchar(refseq))
  for (i in 1:nrow(regions)) {
    col[(regions[i, "start"] + 1):regions[i, "end"]] <- colors[i]
  }
  col <- col[target_start + 1:target_end]
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
#' @param traceQC_input A TraceQC object
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
plot_deletion_hotspot <- function(traceQC_input, count_cutoff = 1) {
  deletions <- traceQC_input$mutation %>%
    filter(type == "deletion") %>%
    group_by(start, length) %>%
    summarise(count = sum(count)) %>%
    mutate(end = start + length) %>%
    ungroup %>%
    mutate(id = 1:n()) %>%
    mutate(log10_count = log10(count))

  circular_chordgram(df = deletions,
                     title = "Deletions",
                     traceQC_input,
                     count_cutoff)
}

#' Display a circos plot that shows overall insertion pattern across the barcodes.
#'
#' @param traceQC_input A TraceQC object
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
  function(traceQC_input, count_cutoff = 1) {
    insertions <- traceQC_input$mutation %>%
      filter(type == "insertion") %>%
      group_by(start, length) %>%
      summarise(count = sum(count)) %>%
      mutate(end = start + length) %>%
      ungroup %>%
      mutate(id = 1:n()) %>%
      mutate(log10_count = log10(count))

    circular_chordgram(df = insertions,
                       title = "Insertions",
                       traceQC_input,
                       count_cutoff)
  }

#' Display a mutation hotspot circos plot.
#'
#' The circos plot shows the frequency of mutation events for each nucleotide.
#'
#' @param traceQC_input A traceQC object
#'
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @return It won't return any specific object.
#'
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_point_mutation_hotspot(example_obj)
#'
plot_point_mutation_hotspot <- function(traceQC_input) {
  mutations <- filter(traceQC_input$mutation, type == "mutation") %>%
    group_by(start, length, mutate_to) %>%
    summarise(count = sum(count)) %>%
    ungroup %>%
    group_by(start) %>%
    mutate(mutate_to = as.character(mutate_to)) %>%
    arrange(count) %>%
    mutate(y = cumsum(count)) %>%
    ungroup
  circular_histogram(mutations, "Mutations", traceQC_input)
}

