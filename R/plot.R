library(readr)
library(ggplot2)

plot_construct <- function(ref) {
  p <- data.frame(text=unlist(strsplit(ref$seq,"")),
                  pos=1:nchar(ref$seq),
                  x=(1:nchar(ref$seq)-1) %% 50,
                  y=-((1:nchar(ref$seq)-1) %/% 50)) %>%
    mutate(region=case_when(pos<ref$barcode_start ~ "adapter",
                  pos>=ref$spacer_start&pos<ref$spacer_end ~ "spacer",
                  pos>=ref$spacer_end&pos<ref$spacer_end+3 ~ "PAM",
                  pos>=ref$barcode_end ~ "adapter"))
    
  ggplot(p) + 
    geom_text(aes(x=x,y=y,label=text,color=region)) +
    coord_fixed(ratio=nchar(ref$seq) %/% 50) +
    theme_void()
}

lorenz_curve <- function(aligned_barcode) {
  p <- aligned_barcode %>%
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
    group_by(barcode) %>%
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
