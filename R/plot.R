library(readr)
library(ggplot2)
library(stringr)
library(DECIPHER)
library(RColorBrewer)

plot_construct <- function(ref) {
  ref <- readLines(ref)
  refseq <- ref[1]
  regions <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions) %>% 
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(start),
           end=strtoi(end)) %>%
    mutate(region=as.character(region))
  
  colors <- brewer.pal(nrow(regions)+1,"Set2")
  p <- data.frame(text=unlist(strsplit(refseq,"")),
                  pos=1:nchar(refseq),
                  x=(1:nchar(refseq)-1) %% 50,
                  y=-((1:nchar(refseq)-1) %/% 50))
  p$region="adapter"
  for (i in 1:nrow(regions)) {
    p$region[(regions[i,"start"]+1):regions[i,"end"]] <- regions[i,"region"]
  }
  
  ggplot(p) + 
    geom_text(aes(x=x,y=y,label=text,color=region)) +
    scale_color_manual(values=colors,breaks=c(regions$region,"adapter")) +
    coord_fixed(ratio=nchar(refseq) %/% 50) +
    theme_void()
}

lorenz_curve <- function(aligned_reads) {
  p <- aligned_reads %>%
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
  p <- mutation_df %>%
    group_by(type,start,length) %>%
    summarise(num_mutation=n()) %>%
    ungroup %>%
    mutate(type=as.factor(type))
  
  p <- ggplot(p) + 
    geom_bar(aes(x=type,fill=type)) +
    labs(x="",y="diversity") +
    theme_classic()
  return(p)
}
