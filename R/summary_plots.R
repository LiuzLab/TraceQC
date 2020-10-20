

#' Visualization of the construct (reference sequence) information.
#'
#' @param traceQC_input an TraceQC object
#'
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#'
#' @return it returns A ggplot2 object that shows the construct information.
#' @export
#'
#' @examples
#' library(TraceQC)
#' data(example_obj)
#' plot_construct(example_obj)
#'
plot_construct <- function(traceQC_input) {
  colors <- brewer.pal(length(unique(traceQC_input$regions$region)) + 1, "Set2")
  p <- data.frame(
    text = unlist(strsplit(traceQC_input$refseq, "")),
    pos = 1:nchar(traceQC_input$refseq),
    x = (1:nchar(traceQC_input$refseq) - 1) %% 50,
    y = -((1:nchar(
      traceQC_input$refseq
    ) - 1) %/% 50)
  )
  p$region <- "adapter"
  for (i in 1:nrow(traceQC_input$regions)) {
    from <- traceQC_input$regions[i, "start"]
    to <- traceQC_input$regions[i, "end"]
    p$region[from:to] <-
      traceQC_input$regions[i, "region"]
  }

  ggplot(p) +
    geom_text(aes_string(
      x = "x",
      y = "y",
      label = "text",
      color = "region"
    )) +
    scale_color_manual(values = colors,
                       breaks = c(unique(traceQC_input$regions$region), "adapter")) +
    coord_fixed(ratio = nchar(traceQC_input$refseq) %/% 50, clip="off") +
    theme_void()
}

#' Drawing a score distribution plot
#'
#' @param traceQC_input A TraceQC object
#'
#' @import ggplot2
#'
#' @return A ggplot2 object that shows alignment score distribution.
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_score_distribution(example_obj)
#'
plot_score_distribution <- function(traceQC_input) {
  ggplot(traceQC_input$aligned_reads) +
    geom_histogram(aes_string(x = "score", y = "..density.."), binwidth = 5) +
    ylab("percentage") +
    theme_classic()
}

#' Drawing Lorenz Curve
#'
#' The Lorenz curve shows an inequality of barcode distribution of the sample.
#'
#' @param traceQC_input A TraceQC object
#'
#' @import ggplot2
#'
#' @return A ggplot2 object that shows Lorenz Curve
#' @export
#'
#' @examples
#' data(example_obj)
#' plot_lorenz_curve(example_obj)
#'
plot_lorenz_curve <- function(traceQC_input) {
  p <- traceQC_input$aligned_reads %>%
    group_by(.data$target_seq) %>%
    summarise(count = n()) %>%
    ungroup %>%
    arrange(desc(count)) %>%
    mutate(x = 1:n(), cs = cumsum(count) / sum(count))

  scale <- log10(p$count[1])
  p <- ggplot(p) +
    geom_line(aes_string(x = "x", y = "y"),
              data = data.frame(x = c(0, nrow(p)), y = c(0, 1)),
              linetype = "dotted") +
    geom_line(aes_string(x = "x", y = "log(count) / scale"), color = "blue") +
    geom_line(aes_string(x = "x", y = "cs"), color = "red") +
    coord_fixed(ratio = nrow(p)) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      sec.axis = sec_axis( ~ . * scale, name = "Log barcode count (blue)")
    ) +
    labs(x = "Ranked barcode", y = "Cumulative sum (red)") +
    theme_classic()
  return(p)
}

#' A barplot to show distribution of the number of mutations per barcode
#'
#' @param traceQC_input A TraceQC object
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A ggplot2 object that shows the barplot
#' @export
#'
#' @examples
#' data(example_obj)
#' num_mutation_histogram(example_obj)
#'
num_mutation_histogram <- function(traceQC_input) {
  p <- traceQC_input$mutation %>%
    group_by(.data$target_seq) %>%
    summarise(num_mutation = n()) %>%
    ungroup %>%
    mutate(num_mutation = case_when(.data$num_mutation >= 10 ~ "10+",
                                    TRUE ~ as.character(.data$num_mutation))) %>%
    mutate(num_mutation = factor(.data$num_mutation, levels = c(as.character(1:9), "10+")))

  p <- ggplot(p) +
    geom_bar(aes_string(x = "num_mutation")) +
    labs(x = "number of mutations per barcode") +
    theme_classic()

  return(p)
}

#' A pie chart that shows a summary of mutation types.
#'
#' @param traceQC_input A TraceQC object
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A ggplot2 object that shows the pie chart
#' @export
#'
#' @examples
#' data(example_obj)
#' mutation_type(example_obj)
#'
mutation_type <- function(traceQC_input) {
  breaks <- c(1, 2, 4, 8, 16)

  df <- traceQC_input$mutation %>%
    filter(.data$type != "unmutated") %>%
    group_by(.data$type, .data$start, .data$length, .data$mutate_to) %>%
    summarise(count = n()) %>%
    ungroup %>%
    mutate(length_category = findInterval(.data$length, breaks)) %>%
    group_by(.data$type, .data$length_category) %>%
    summarise(count = n()) %>%
    ungroup %>%
    arrange(.data$type, .data$length_category) %>%
    mutate(ymax = cumsum(.data$count) / sum(.data$count),
           labels = breaks[2:length(breaks)][length_category])

  plotting_df <- df %>%  group_by(.data$type) %>%
    summarise(ymax = max(.data$ymax), count = sum(.data$count)) %>%
    ungroup

  plotting_df$ymin <- c(0, plotting_df$ymax[1:2])
  plotting_df$labelPosition = (plotting_df$ymax + plotting_df$ymin) / 2
  plotting_df$label <-
    paste0(plotting_df$type, "\n value: ", plotting_df$count)

  p <-
    ggplot(plotting_df,
           aes_string(
             ymax = "ymax",
             ymin = "ymin",
             xmax = "4",
             xmin = "3",
             fill = "type"
           )) +
    geom_rect() +
    geom_label(x = 3.5,
               aes_string(y = "labelPosition", label = "label"),
               size = 3) +
    scale_fill_brewer(palette = 4) +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) +
    scale_y_continuous(
      breaks = df %>% filter(.data$type != "mutation") %>% pull(.data$ymax),
      labels = df %>% filter(.data$type != "mutation") %>% pull(.data$labels)
    ) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12)
    )
  return(p)
}
