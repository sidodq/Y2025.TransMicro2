#' @title Abundance table : taxonomy with multi-affiliation and abundance file without chimeras
#' @param multi_affi_table multi-affiliation table cleaned
#' @param abundance_table abundance table cleaned
#' @param name name will be used to import files and write names before extension of the different output tables (for function find_chimeras : _cleaned_abundance.tsv, _chimeras_list.tsv, _chimeras_summary.txt)
#'
#' @examples
#' at_MS_16S <- Combine_ab_multi_affi_tables_remove_chim(multi_affi_table = taxo_16S_Multiaffiliation, name = "MicrobiotaStructure_16S", abundance_table = ab_MS_16S)
#'
#' @description
#' - Merge abundance table taxonomy (y=blast_taxonomy) and multi-affiliation table taxonomy (x=multi_affi_table)
#' - Create a column blast_taxonomy keeping the multi-affiliation taxonomy but if NA then take the abundance taxonomy
#' - Identify and remove chimeras from the abundance table (with function Remove_chimeras)
#' @return
#' An abundance table with a clean taxonomy (with accounting for a threshold for multi-affiliation) and without chimeras
#'
#' @export
Combine_ab_multi_affi_tables_remove_chim <- function(abundance_table, multi_affi_table, name)

{
  # Merge abundance table taxonomy (x=blast_taxonomy) and multi-affiliation table taxonomy (y=multi_affi_table)
  my.result <- merge(multi_affi_table, abundance_table, by="observation_name", all.y=TRUE)
  abundance_table$blast_taxonomy
  # Create a column blast_taxonomy keeping the multi-affiliation taxonomy (x) but if NA then take the abundance taxonomy (y)
  my.result$blast_taxonomy <- ifelse(!is.na(my.result$blast_taxonomy.x), my.result$blast_taxonomy.x, my.result$blast_taxonomy.y)
  # Deleting columns
  my.result$blast_taxonomy.x<-NULL
  my.result$blast_taxonomy.y<-NULL
  # Reordering columns
  data.table::setcolorder(my.result, c(2,3,4,5,6,7,8,9,10,1,11:(ncol(my.result))))
  fordada <- as.data.frame(my.result)
  fordada <- fordada[order(as.integer(gsub("Cluster_", "", fordada$observation_name))),]

  # Remove chimeras and export
  Remove_chimeras(fordada=fordada, name=name)
}
