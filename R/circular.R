library(tidyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)

#' Title
#'
#' @param df
#' @param title
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
circular_chordgram <- function(df,title,ref) {
  ref <- readLines(ref)
  regions <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions) %>%
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(start),
           end=strtoi(end)) %>%
    mutate(region=as.character(region))
  target_start <- filter(regions,region=="target") %>% pull(start)
  target_end <- filter(regions,region=="target") %>% pull(end)
  refseq <- substr(ref[1],start=target_start+1,stop=target_end)

  col_fun = colorRamp2(c(floor(min(df$log10_count)),
                         ceiling(max(df$log10_count))),c("yellow", "red"))
  df$color <- col_fun(df$log10_count)
  l <- nchar(refseq) + 1
  circos.initialize(factors = factor(1), xlim = c(1, l))
  circos.track(track.height=0.1,ylim=c(0,1))
  circos.axis(major.at=seq(0,l,10),labels=seq(0,l,10))
  for (i in 1:nrow(df)) {
  circos.link(1,df$start[i],1,df$end[i],h.ratio=0.9,
              lwd=0.2*df$log10_count[i],col=df$color[i])
  }

  colors <- brewer.pal(nrow(regions)+1,"Set2")
  col = rep("black",nchar(ref[1]))
  for (i in 1:nrow(regions)) {
    col[(regions[i,"start"]+1):regions[i,"end"]] <- colors[i]
  }
  col <- col[target_start+1:target_end]

  circos.text(1:l,0.5,strsplit(refseq,split="") %>% unlist(),col=col,cex=1)
  title(title)
  circos.clear()
  lgd <- Legend(at=seq(floor(min(df$log10_count)),ceiling(max(df$log10_count)),length.out=5),
                col_fun=col_fun,title="log 10 count")
  draw(lgd,x = unit(0.15, "npc"), y = unit(0.15, "npc"))
}


#' Title
#'
#' @param df
#' @param ref
#'
#' @return
#' @export
#'
#' @examples
circular_histogram <- function(df,ref) {
  ref <- readLines(ref)
  regions <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions) %>%
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(start),
           end=strtoi(end)) %>%
    mutate(region=as.character(region))
  target_start <- filter(regions,region=="target") %>% pull(start)
  target_end <- filter(regions,region=="target") %>% pull(end)
  refseq <- substr(ref[1],start=target_start+1,stop=target_end)

  l <- nchar(refseq) + 1
  circos.initialize(factors = factor(1), xlim = c(1, l))
  circos.track(track.height=0.2,ylim=c(0,300),bg.border=NA)
  circos.track(track.height=0.1,ylim=c(0,1))
  circos.axis(h="bottom",labels.facing="reverse.clockwise",direction="inside",
              major.at=seq(0,l,10),labels=seq(0,l,10))
  colors <- brewer.pal(nrow(regions)+1,"Set2")
  col = rep("black",nchar(ref[1]))
  for (i in 1:nrow(regions)) {
    col[(regions[i,"start"]+1):regions[i,"end"]] <- colors[i]
  }
  col <- col[target_start+1:target_end]
  circos.text(1:l,0.5,strsplit(refseq,split="") %>% unlist(),col=col,cex=1)

  colors <- c("red","grey","blue","green")
  names(colors) <- c("A","C","G","T")
  for (i in 1:nrow(df)) {
    circos.rect(df$start[i]-0.4,
                df$y[i]-df$count[i]+1,
                df$start[i]+0.4,
                df$y[i],
                col=colors[df$mutate_to[i]],
                border=NA,track.index=1)
  }
  circos.clear()}
