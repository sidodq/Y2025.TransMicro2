#' @title Charge 3 xlsx to form a physeq object
#'
#' @description
#' Function to create a physeq object from three xlsx docs (data, otu & tax tables)
#'
#' @param ps_name ex. "ps_rare_ITS_D1"
#' @param folder_path here::here() path to the FOLDER where to take the files
#'
#' @export
Charge_xl_to_ps <- function(ps_name, folder_path) {

  #Charge
  otu_path <- file.path(folder_path, paste0(ps_name, "_otu_table.xlsx"))
  otu <- read.xlsx(xlsxFile=otu_path, rowNames = TRUE)

  rownames(otu) <- otu$observation_name
  otu$observation_name <- NULL
  otu[] <- lapply(otu, as.numeric)
  otu_mat <- t(as.matrix(otu))

  OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)



  tax_path <- file.path(folder_path, paste0(ps_name, "_tax_table.xlsx"))
  tax <- read.xlsx(xlsxFile=tax_path, rowNames = TRUE)

  rownames(tax) <- tax$observation_name
  tax$observation_name <- NULL

  TAX <- phyloseq::tax_table(as.matrix(tax))



  mtd_path <- file.path(folder_path, paste0(ps_name, "_mtd_table.xlsx"))
  mtd <- read.xlsx(xlsxFile=mtd_path, rowNames = TRUE)

  MTD <- phyloseq::sample_data(as.data.frame(mtd))

  # Build phyloseq object
  ps <- phyloseq::phyloseq(OTU, TAX, MTD)

  return(ps)
}
