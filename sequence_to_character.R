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

seq_to_character <- function(aligned_reads,out_file,cores=5) {
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
                        target_seqs,read_counts,SIMPLIFY=FALSE,mc.cores=cores) %>%
    bind_rows()

  unmutated <- filter(aligned_reads,target_seq==target_ref)
  mutation_df <- rbind(data.frame(target_seq=unmutated$target_seq,type="unmutated",start=0,
                                  length=0,mutate_to="-",count=unmutated$count),mutation_df)

  write_tsv(mutation_df,out_file)
  return(mutation_df)}
