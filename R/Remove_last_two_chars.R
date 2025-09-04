#' @title Delete the last two characters of a word
#'
#' @description
#' Function to delete the last two characters of a word
#' - For example : Delete the last two characters of the sample names that distinguish the replicas, so that they can then be merged
#'
#' @param text character string
#'
#' @examples
#' data_MS_ITS$Samples <- sapply(data_MS_ITS$Barcode, Remove_last_two_chars)
#'
#' @return
#' An character string with two characters less in each name.
#'


#' @export
Remove_last_two_chars <- function(text) {
  gsub("_(A|B)$", "", text)
}
