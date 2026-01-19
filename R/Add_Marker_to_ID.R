#' @title Add the marker at the end of the sample ID to distinguish 16S replicates from ITS
#'
#' @description
#'
#'
#' @param medata dataset
#' @param samples list of samples names
#'
#' @examples
#' samples <- c("CT_INDEX_3_rn2", "CT_INDEX_5_rn1")
#'
#' @return
#'
#' @export
add_marker_to_id <- function(metadata, samples) {
  metadata %>%
    mutate(
      ID_pool = if_else(
        ID_pool %in% samples,
        paste0(ID_pool, "_", Marker),
        ID_pool
      )
    )
}
