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
Remove_last_two_chars_run <- function(text) {

  # Step 1: if ends with _rnX_Y → swap to _Y_rnX
  text <- gsub("_(rn[12])_([123])$", "_\\2_\\1", text)

  # Step 2: remove trailing _rn1 or _rn2
  text <- gsub("_(rn[12])$", "", text)

  return(text)
}
