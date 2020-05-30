library(readr)
library(ggplot2)
library(stringr)
library(DECIPHER)
library(RColorBrewer)

plot_construct <- function(traceQC_input) {
  colors <- brewer.pal(nrow(traceQC_input$regions)+1,"Set2")
  p <- data.frame(text=unlist(strsplit(traceQC_input$refseq,"")),
                  pos=1:nchar(traceQC_input$refseq),
                  x=(1:nchar(traceQC_input$refseq)-1) %% 50,
                  y=-((1:nchar(traceQC_input$refseq)-1) %/% 50))
  p$region="adapter"
  for (i in 1:nrow(traceQC_input$regions)) {
    p$region[(traceQC_input$regions[i,"start"]+1):traceQC_input$regions[i,"end"]] <- traceQC_input$regions[i,"region"]
  }

  ggplot(p) +
    geom_text(aes(x=x,y=y,label=text,color=region)) +
    scale_color_manual(values=colors,breaks=c(traceQC_input$regions$region,"adapter")) +
    coord_fixed(ratio=nchar(traceQC_input$refseq) %/% 50) +
    theme_void()
}

plot_score_distribution <- function(traceQC_input) {
  ggplot(traceQC_input$aligned_reads) +
    geom_histogram(aes(x=score,y=..density..),binwidth=5) +
    ylab("percentage") +
    theme_classic()
}

lorenz_curve <- function(traceQC_input) {
  p <- traceQC_input$aligned_reads %>%
    group_by(target_seq) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(desc(count)) %>%
    mutate(x=1:n(),cs=cumsum(count)/sum(count))

  scale <- log(p$count[1])
  p <- ggplot(p) +
    geom_line(aes(x=x,y=y),data=data.frame(x=c(0,nrow(p)),y=c(0,1)),linetype="dotted") +
    geom_line(aes(x=x,y=log(count)/scale),color="blue") +
    geom_line(aes(x=x,y=cs),color="red") +
    coord_fixed(ratio=nrow(p)) +
    scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1),
                       sec.axis = sec_axis(~.*scale,name="log barcode count")) +
    labs(x="ranked barcode",y="cumulative sum") +
    theme_classic()
  return(p)
}

num_mutation_histogram <- function(mutation_df) {
  p <- mutation_df %>%
    group_by(target_seq) %>%
    summarise(num_mutation=n()) %>%
    ungroup %>%
    mutate(num_mutation=case_when(num_mutation>=10 ~ "10+",
                  TRUE ~ as.character(num_mutation))) %>%
    mutate(num_mutation=factor(num_mutation,levels= c(as.character(1:9),"10+")))

  p <- ggplot(p) +
    geom_bar(aes(x=num_mutation)) +
    labs(x="number of mutations per barcode") +
    theme_classic()

  return(p)
}

mutation_type <- function(mutation_df) {
  breaks <- c(1,2,4,8,16)

  df <- mutation_df %>%
    filter(type!="unmutated") %>%
    group_by(type,start,length,mutate_to) %>%
    summarise(count=n()) %>%
    ungroup %>%
    mutate(length_category=findInterval(length,breaks)) %>%
    group_by(type,length_category) %>%
    summarise(count=n()) %>%
    ungroup %>%
    arrange(type,length_category) %>%
    mutate(ymax=cumsum(count)/sum(count))

  plotting_df <- group_by(df,type) %>%
    summarise(ymax=max(ymax),count=sum(count)) %>%
    ungroup

  plotting_df$ymin <- c(0,plotting_df$ymax[1:2])
  plotting_df$labelPosition = (plotting_df$ymax+plotting_df$ymin)/2
  plotting_df$label <- paste0(plotting_df$type, "\n value: ", plotting_df$count)

  p <- ggplot(plotting_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=type)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
    scale_fill_brewer(palette=4) +
    coord_polar(theta="y") +
    xlim(c(2, 4)) +
    scale_y_continuous(breaks=filter(df,type!="mutation",length_category!=5) %>% pull(ymax),
                       labels=c(breaks[2:length(breaks)],breaks[2:length(breaks)])) +
    theme(legend.position = "none",
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color="grey"),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.title = element_blank(),
          axis.text.x = element_text(face="bold",size=12))
  return(p)
}
