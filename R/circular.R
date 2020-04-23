library(tidyr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)

circular_chordgram <- function(df,title,ref){
  ref_seq <- substr(ref$seq,start=ref$barcode_start,stop=ref$barcode_end-1)
  col_fun = colorRamp2(c(floor(min(df$log10_count)),
                         ceiling(max(df$log10_count))),c("yellow", "red"))
  df$color <- col_fun(df$log10_count)
  l <- nchar(ref_seq) + 1
  circos.initialize(factors = factor(1), xlim = c(1, l))
  circos.track(track.height=0.1,ylim=c(0,1))
  circos.axis(major.at=seq(0,l,10),labels=seq(0,l,10))
  for (i in 1:nrow(df)) {
  circos.link(1,df$start[i],1,df$end[i],h.ratio=0.9,
              lwd=0.2*df$log10_count[i],col=df$color[i])
  }
  col <- c(rep("black",ref$spacer_start-ref$barcode_start), 
           rep("blue",ref$spacer_end-ref$spacer_start), 
           rep("green",3), 
           rep("black",ref$barcode_end-ref$spacer_end-3))
  circos.text(1:l,0.5,strsplit(ref_seq,split="") %>% unlist(),col=col,cex=1)
  title(title)
  circos.clear()
  lgd <- Legend(at=seq(floor(min(df$log10_count)),ceiling(max(df$log10_count)),length.out=5),
                col_fun=col_fun,title="log 10 count")
  draw(lgd,x = unit(0.15, "npc"), y = unit(0.15, "npc"))
}

circular_histogram <- function(df,ref) {
  ref_seq <- substr(ref$seq,start=ref$barcode_start,stop=ref$barcode_end-1)
  l <- nchar(ref_seq) + 1
  scale <- df %>% 
    group_by(start) %>% 
    summarise(count=sum(count)) %>%
    ungroup
  regions <- c(rep("adapter",ref$spacer_start-ref$barcode_start), 
           rep("spacer",ref$spacer_end-ref$spacer_start), 
           rep("PAM",3), 
           rep("adapter",ref$barcode_end-ref$spacer_end-3))
  text_df <- data.frame(x=1:(l-1),y=-0.2*max(scale$count),
                label=strsplit(ref_seq,split="") %>% unlist(),regions=regions)
  
  ggplot(df) +
    geom_bar(aes(x=start,y=count,fill=mutate_to),stat="identity") +
    geom_text(aes(x=x,y=y,label=label,color=regions),data=text_df) +
    scale_color_manual(values=c("adapter"="black","spacer"="blue","PAM"="green")) +
    coord_polar(start = pi/2) +
    theme(panel.grid = element_blank(),panel.background=element_blank(),
          axis.title=element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_x_continuous(breaks=seq(0,l,10),labels=seq(0,l,10),lim=c(0,l+5)) +
    scale_y_continuous(lim=c(-5*max(scale$count),1.2*max(scale$count)))}



