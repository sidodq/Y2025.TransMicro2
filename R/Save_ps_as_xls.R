#' @title Save physeq object as 3 readable files
#'
#' @description
#' Function to export physeq object as readable excel docs (data, otu & tax tables)
#'
#' @param ps phyloseq object
#' @param file_path here::here() path to the FOLDER where to save the files
#'
#'
#' @export
Save_ps_as_xls <- function(ps, folder_path) {

  ps_name <- deparse(substitute(ps))

  otu_df <- as.data.frame(otu_table(ps))
  otu_df$observation_name <- rownames(otu_df)
  otu_df<- otu_df %>% dplyr::relocate(observation_name)

  otu_file <- file.path(folder_path, paste0(ps_name, "_otu_table.xlsx"))

  write.xlsx(otu_df, file = otu_file, rowNames = TRUE)



  tax_df <- as.data.frame(tax_table(ps))
  tax_df$observation_name <- rownames(tax_df)
  tax_df <- tax_df %>% dplyr::relocate(observation_name)

  tax_file <- file.path(folder_path, paste0(ps_name, "_tax_table.xlsx"))

  write.xlsx(tax_df, file= tax_file, rowNames = TRUE)



  metadata_df <- data.frame(sample_data(ps), check.names = FALSE)

  mtd_file <- file.path(folder_path, paste0(ps_name, "_mtd_table.xlsx"))

  write.xlsx(tax_df, file= mtd_file, rowNames = TRUE)
}
