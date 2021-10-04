#' Get read count per UMI.
#'
#' @param df aligned reads
#'
#' @return A data frame of UMI
#' @export
#'
get_read_count_per_UB <- function(df) {
  group_by(df,UB,CB) %>%
    summarise(read_count_per_UMI=n()) %>%
    ungroup
}


#' Get UMI count per cell.
#'
#' @param df aligned reads
#'
#' @return A data frame of Cells.
#' @export
#'
get_UMI_count_per_CB <- function(df) {
  group_by(df,CB,UB) %>%
    summarise() %>%
    ungroup %>%
    group_by(CB) %>%
    summarise(UMI_per_CB=n()) %>%
    ungroup
}

#' Filter mutations based on read count per UMI
#'
#' @param data A data frame.
#' @param freq_threshold threshold of mutation frequency.
#' @param include_max include the mutations with maximum read count.
#'
#' @return A filtered data frame.
#' @export
#'
filter_mutations <- function(data,include_max=TRUE,freq_threshold) {
  df <- filter(data,count>freq_threshold*read_count_per_UMI)
  if (!include_max) {return(df)} else {
    if (nrow(df)>0) {
      return(df)
    }
    max_count <- max(data$count)
    return(filter(data,count==max_count))
  }
}

#' Filter mutations based on read count per UMI
#'
#' @param data A data frame.
#' @param freq_threshold threshold of mutation frequency.
#' @param include_max include the mutations with maximum read count.
#'
#' @importFrom tidyr nest 
#' @importFrom tidyr unnest
#' @importFrom purrr map
#' 
#' @return A filtered data frame.
#' @export
#'
filter_mutations_per_UMI <- function(data,include_max=TRUE,freq_threshold) {
  df <- group_by(data,CB,UB,type,start,length,mutate_to) %>%
    summarise(count=sum(count),read_count_per_UMI=mean(read_count_per_UMI)) %>%
    ungroup
  
  group_by(df,CB,UB) %>%
    nest() %>%
    mutate(mutations=map(data,filter_mutations,include_max,freq_threshold)) %>%
    unnest(mutations) %>%
    select(-data) %>%
    ungroup}


#' Filter mutations based on read count per cell
#'
#' @param data A data frame.
#' @param freq_threshold threshold of mutation frequency.
#'
#' @importFrom dplyr filter 
#' @return A filtered data frame.
#' @export
#'
filter_mutations_per_cell <- function(data,freq_threshold) {
  group_by(data,CB,type,start,length,mutate_to) %>%
    summarise(UMI_count=n(),UMI_per_CB=mean(UMI_per_CB)) %>%
    ungroup %>%
    filter(UMI_count >= freq_threshold*UMI_per_CB)}
