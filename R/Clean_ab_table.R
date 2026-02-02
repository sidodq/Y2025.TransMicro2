#' @title Clean abundance table
#' @param abundance_table Abundance table obtained from FROGS
#' @param remove_col Remove columns containing this specific word
#' @description
#'- Eliminate specific columns/samples (example : "Undetermined")
#'- Cleaning up the taxonomy using function Clean_taxo
#'
#' @return
#' An abundance table with a clean taxonomy (but without accounting for a threshold for multi-affiliation) and with the desired samples
#'
#' @examples
#' MicrobiotaStructure_16S_Abundance <- Clean_taxo_at(abundance_table = at, remove_col = NULL)
#'
#' @export
Clean_ab_table <- function(abundance_table, remove_col)

{
  # Delete columns containing the word indicated in remove_col and display the names of the remaining columns
  if(!is.null(remove_col)){
    abundance_table <- abundance_table %>%
      dplyr::select(-contains(remove_col))
    print("Remaining colnames")
    print(colnames(abundance_table))
  }

  # Create bidouille column (a copy of the blast_taxonomy column)
  ab <- data.table::data.table(bidouille=abundance_table$`blast_taxonomy`, abundance_table)

  # Clean the taxonomy from blast in ab (bidouille column)
  my_ab_result <- Clean_taxo(table = ab)

  # Reconstruction of the clean full taxonomy from each taxonomic field : assemble taxonomic ranks using ;
  my_ab_result$blast_taxonomy <- apply(my_ab_result, 1, function(r) {
    paste(r[1:7], sep=";", collapse=";")
  })

  # Delete columns tx1, tx2, tx3, tx4, tx5, tx6, tx7
  my_ab_result <- subset(my_ab_result, select = -c(tx1, tx2, tx3, tx4, tx5, tx6, tx7))

  return(my_ab_result)
}

