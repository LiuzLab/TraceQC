seq_to_character <- function(aligned_barcode) {
  all_insertions <- str_locate_all(aligned_barcode$ref_aligned,"-+")
  all_deletions <- str_locate_all(aligned_barcode$seq_aligned,"-+")
  all_mutations <- mapply(function(x,y) {which(x!=y & x!="-" & y!="-" & y!="N")},
                        strsplit(aligned_barcode$ref_aligned,""), strsplit(aligned_barcode$seq_aligned,""))
  
  mutation_df <- data.frame()
  for (i in 1:nrow(aligned_barcode)) {
  insertions <- all_insertions[[i]] %>%
    as.data.frame() %>%
    mutate(barcode=aligned_barcode$seq_aligned[i]) %>%
    mutate(length=end-start+1,type="insertion",mutate_to=substr(barcode,start=start,stop=end))
  deletions <- all_deletions[[i]] %>%
    as.data.frame() %>%
    mutate(barcode=aligned_barcode$seq_aligned[i]) %>%
    mutate(length=end-start+1,type="deletion",mutate_to="-")
  mutations <- data.frame(start=all_mutations[[i]],end=all_mutations[[i]]) %>%
    mutate(barcode=aligned_barcode$seq_aligned[i]) %>%
    mutate(length=1,type="mutation",mutate_to=substr(barcode,start=start,stop=end))
  
  mutation_df <- rbind(mutation_df,rbind(insertions,deletions,mutations) %>%
                         arrange(start) %>%
                         mutate(tmp=lag(as.integer(type=="insertion")*length)) %>%
                         mutate(tmp=replace_na(tmp,0)) %>%
                         mutate(align=cumsum(tmp)) %>%
                         mutate(start=start-align,count=aligned_barcode$count[i]) %>%
                         select(count,barcode,type,start,length,mutate_to))}
  write_tsv(mutation_df,"./mutational_character.tsv")}
