#' Visualization of the construct (reference sequence) information.
#'
#' @param ref an reference object, output of `parse_ref_file`
#' @param chr_per_row number of charachters per row.
#' @param chr_size the font size of characters.
#'
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#'
#' @return it returns A ggplot2 object that shows the construct reference sequence.
#' @export
#' @examples 
#' 
#' ref_file <- system.file("extdata/test_data/ref","ref_carlin.txt",package="TraceQC")
#' ref <- parse_ref_file(ref_file)
#' plot_construct(ref,chr_per_row=50,chr_size=5)
#'
plot_construct <- function(ref,chr_per_row=50,chr_size=10) {
  colors <- brewer.pal(length(unique(ref$regions$region)) + 1, "Set2")
  p <- data.frame(
    text = unlist(strsplit(ref$refseq, "")),
    pos = 1:nchar(ref$refseq),
    x = (1:nchar(ref$refseq) - 1) %% chr_per_row,
    y = -((1:nchar(
      ref$refseq
    ) - 1) %/% chr_per_row)
  )
  p$region <- "adapter"
  for (i in 1:nrow(ref$regions)) {
    from <- ref$regions[i, "start"]
    to <- ref$regions[i, "end"]
    p$region[from:to] <-
      ref$regions[i, "region"]
  }

  ggplot(p) +
    geom_text(aes_string(
      x = "x",
      y = "y",
      label = "text",
      color = "region"
    ),size=chr_size) +
    scale_color_manual(values = colors,
                       breaks = c(unique(ref$regions$region), "adapter")) +
    coord_fixed(ratio = nchar(ref$refseq) %/% chr_per_row, clip="off") +
    theme_void()
}

#' Visualization of alignment permutation.
#'
#' @param ref an data frame of permutation sequence, output of `sequence_permutation`
#'
#' @import ggplot2
#'
#' @return it returns A ggplot2 object that shows the permutation.
#' @export
#'
plot_alignment_permutation <-  function(alignment_permutation) {
    model <- loess(score~permutate_percent,data=alignment_permutation)
    ggplot(alignment_permutation,aes(x=permutate_percent,y=score)) +
      geom_point() +
      geom_smooth(formula=y~x,method="loess",sd=TRUE) +
      theme_classic()}

#' Drawing a score distribution plot
#'
#' @param aligned_reads A aligned_reads dataframe.
#'
#' @import ggplot2
#'
#' @return A ggplot2 object that shows alignment score distribution.
#' @export
#'
#' @examples
#' plot_score_distribution(aligned_reads)
#'
plot_score_distribution <- function(aligned_reads) {
  ggplot(aligned_reads) +
    geom_histogram(aes_string(x = "score", y = "..density.."), binwidth = 5) +
    ylab("percentage") +
    theme_classic()
}

#' Drawing Lorenz Curve
#'
#' The Lorenz curve shows an inequality of barcode distribution of the sample.
#'
#' @param aligned_reads A aligned_reads dataframe.
#'
#' @import ggplot2
#'
#' @return A ggplot2 object that shows Lorenz Curve
#' @export
#'
#' @examples
#' plot_lorenz_curve(aligned_reads)
#'
plot_lorenz_curve <- function(aligned_reads) {
  p <- aligned_reads %>%
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
#' @param mutations A mutation dataframe
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A ggplot2 object that shows the barplot
#' @export
#'
#' @examples
#' mutation <- read_tsv(system.file("extdata/test_data","mutation.txt",package="TraceQC"))
#' num_mutation_histogram(mutation)
#'
num_mutation_histogram <- function(mutations) {
  p <- mutations %>%
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
#' @param mutations A mutation dataframe
#'
#' @import ggplot2
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @return A ggplot2 object that shows the pie chart
#' @export
#'
#' @examples
#' mutation <- read_tsv(system.file("extdata/test_data","mutation.txt",package="TraceQC"))
#' mutation_type_donut(mutation)
#'
mutation_type_donut <- function(mutations) {
  plotting_df <- mutations %>%
    filter(.data$type != "unmutated") %>%
    group_by(.data$type, .data$start, .data$length, .data$mutate_to) %>%
    summarise() %>%
    ungroup %>%
    group_by(.data$type) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(.data$type) %>%
    mutate(ymax = cumsum(.data$count) / sum(.data$count))
  
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
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.background = element_rect(fill = "transparent", colour = NA),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
    )
  return(p)}
