create_input_object <- function(aligned_reads_file,ref_file) {
  aligned_reads <- read_tsv(aligned_reads_file)

  ref <- readLines(ref_file)
  refseq <- ref[1]
  regions <- strsplit(ref[2:length(ref)],split=" ")
  regions <- do.call(rbind,regions) %>%
    as.data.frame() %>%
    setNames(c("region","start","end")) %>%
    mutate(start=strtoi(start),
           end=strtoi(end)) %>%
    mutate(region=as.character(region))
  return (list(aligned_reads=aligned_reads,
               refseq=refseq,
              regions=regions))
}

library(reticulate)
sequence_alignment <- function(input_file,ref_file,output_file="./output/aligned_reads.txt",
                      match=2,mismatch=-2,gapopen=-6,gapextension=-0.1) {
  args <- list("input"=input_file,
               "reference"=ref_file,
               "output"=output_file,
               "match"=match,"mismatch"=mismatch,"gapopen"=gapopen,
               "gapextension"=gapextension)
  source_python(system.file("py", "alignment.py", package="traceQC"))
  alignment(args)
}

#' Title
#'
#' @param insertions
#' @param deletions
#' @param mutations
#' @param target_seq
#' @param read_count
#'
#' @return
#' @export
#'
#' @examples
find_position <- function(insertions,deletions,mutations,target_seq,read_count) {
  insertions <- insertions %>%
    as.data.frame() %>%
    mutate(target_seq=target_seq) %>%
    mutate(length=end-start+1,type="insertion",mutate_to=substr(target_seq,start=start,stop=end))
  deletions <- deletions %>%
    as.data.frame() %>%
    mutate(target_seq=target_seq) %>%
    mutate(length=end-start+1,type="deletion",mutate_to="-")
  mutations <- data.frame(start=mutations,end=mutations) %>%
    mutate(target_seq=target_seq) %>%
    mutate(length=1,type="mutation",mutate_to=substr(target_seq,start=start,stop=end))
  return (rbind(insertions,deletions,mutations) %>%
            arrange(start) %>%
            mutate(tmp=lag(as.integer(type=="insertion")*length)) %>%
            mutate(tmp=replace_na(tmp,0)) %>%
            mutate(align=cumsum(tmp)) %>%
            mutate(start=start-align,
                   target_seq=target_seq,
                   count=read_count) %>%
            select(target_seq,type,start,length,mutate_to,count))
}

#' Title
#'
#' @param aligned_reads
#' @param out_file
#'
#' @return
#' @export
#'
#' @examples
seq_to_character <- function(traceQC_input,out_file) {
  aligned_reads <- traceQC_input$aligned_reads
  col_names <- names(aligned_reads)
  col_names <- col_names[!(col_names %in% c("name","seq","ref","score"))]

  aligned_reads <- aligned_reads %>%
    group_by_at(col_names) %>%
    summarise(count=n()) %>%
    ungroup

  all_insertions <- str_locate_all(aligned_reads$target_ref,"-+")
  all_deletions <- str_locate_all(aligned_reads$target_seq,"-+")
  all_mutations <- mapply(function(x,y) {which(x!=y & x!="-" & y!="-" & y!="N")},
                          strsplit(aligned_reads$target_ref,""), strsplit(aligned_reads$target_seq,""))
  read_counts <- aligned_reads$count
  target_seqs <- aligned_reads$target_seq

  mutation_df <- mcmapply(find_position,all_insertions,all_deletions,all_mutations,
                        target_seqs,read_counts,SIMPLIFY=FALSE,mc.cores=5) %>%
    bind_rows()

  unmutated <- filter(aligned_reads,target_seq==target_ref)
  mutation_df <- rbind(data.frame(target_seq=unmutated$target_seq,type="unmutated",start=0,
                                  length=0,mutate_to="-",count=unmutated$count),mutation_df)

  write_tsv(mutation_df,out_file)
  return(mutation_df)}


build_character_table <- function(df) {
  seq_id <- group_by(df,target_seq,count) %>%
    summarise() %>%
    ungroup %>%
    mutate(sequenceID=1:n())

  unique_events <- df %>%
    group_by(type,start,length,mutate_to) %>%
    summarise() %>%
    ungroup %>%
    mutate(mutationID=1:n())

  nr <- nrow(seq_id)
  nc <- nrow(unique_events)

  df <- left_join(df,seq_id) %>%
    left_join(unique_events) %>%
    mutate(coord=nr*(mutationID-1)+sequenceID)

  character_table <- matrix(data=0,nrow=nr,ncol=nc)
  character_table[df$coord] <- 1
  rownames(character_table) <- 1:nr
  return(list(character_table=character_table,
              seq_id=seq_id))}


