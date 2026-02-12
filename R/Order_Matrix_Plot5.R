#' @title Create heatmap with LBM clusters
#'
#' @param SBMobject BipartiteSBM object
#' @param BF_mat binary matrix
#' @param mtd the corresponding metadata (samples info)
#' @param varibale_name name of the variable for colouring the squares
#' @param modality_string modalities of the variable to colour the squares displayed in the legend
#' @param color_string character string with colours for each mode of the variable
#' @param legend_name explicit legend title
#'
#' @description
#' Creation of a heatmap with samples in y and OTUs in x. A white square represents the absence of the OTU in the sample. A coloured square represents the presence of the OTU in the sample. Samples and OTUs are grouped according to the membership provided by the Latent Block Model. The colour of the squares is defined according to an environmental variable that can be chosen.
#'
#' @return
#' It's a "ggplot" class object.
#'
#' @examples
#' LBM_16S_MPC3_Fruit <- Order_Matrix_Plot(SBMobject = LBM_bernoulli_MPC3_16S, BF_mat = BF_MPC3_16S, varibale_name = "Fruit", modality_string = c("Blackberry", "Cherry", "Strawberry"), color_string = Color_Fruit)

#'
#' @export
Order_Matrix_Plot5 <- function(SBMobject, BF_mat, mtd, modality_string, variable_name, color_string,legend_name) {
  ## Use presence/absence table with OTU as rows and samples as columns
  df = as.data.frame(BF_mat)
  df$taxon = rownames(df)
  ## Transform to long format
  df_long <-
    tidyr::pivot_longer(df,
                        cols = names(df)[-(dim(BF_mat)[[2]]+1)],
                        names_to = "Sample",
                        values_to = "PresenceAbsence")
  ## Add information from metadata
  df_long <- dplyr::left_join(df_long, mtd, by = "Sample")

  ## Number of blocks for rows (samples) and columns (OTUs)
  Qrow <- SBMobject$nbBlocks["row"][[1]]
  Qcol <- SBMobject$nbBlocks["col"][[1]]

  ## Assign OTU cluster (from 1 to Qcol) to each OTU
  df_tax_cluster = data.frame(taxon = rownames(BF_mat),
                              tax_cluster = SBMobject$memberships$col)
  df_long <- dplyr::left_join(df_long, df_tax_cluster, by = "taxon")

  ## Order dataset per tax_cluster
  df_long$tax_cluster <- factor(df_long$tax_cluster, levels = c(1,3,5,2,6,4,7,8))

  df_long <- df_long |>
    dplyr::arrange(tax_cluster)

  ## Order OTU names by tax_cluster
  correct_taxon_order <- unique(df_long$taxon)

  ## Assign sample cluster (from 1 to Qrow) to each sample
  df_spl_cluster = data.frame(Sample = colnames(BF_mat),
                              spl_cluster = SBMobject$memberships$row)

  df_long <- dplyr::left_join(df_long, df_spl_cluster, by = "Sample")

  ## Order dataset per sample_cluster
  df_long$spl_cluster <- factor(df_long$spl_cluster, levels = c(1, 4, 2, 6, 5,3))

  df_long <- df_long |>
    dplyr::arrange(spl_cluster)
  ## Order sample names by sample_cluster
  correct_sample_order <- unique(df_long$Sample)

  ## Convert column of interest as factor
  df_long[[variable_name]] <- factor(df_long[[variable_name]])

  ## Ass color to each level of the factor of interest
  df_long$color_bin <- factor(ifelse(df_long$PresenceAbsence==0, 0, paste(df_long$PresenceAbsence, df_long[[variable_name]], sep = "_")))
  # levels(df_long$color_bin) <- c("", modality_string)
  df_long$color_bin <- ifelse(df_long$color_bin == "0", "", sub("^..", "", df_long$color_bin))
  df_long$color_bin <- factor(df_long$color_bin)

  heatmap_color_bin <- c("white", color_string)

  df_long$PresenceAbsence <- as.factor(df_long$PresenceAbsence)

  # Compute block sizes AFTER reordering
  block_sizes <- df_long |>
    dplyr::distinct(Sample, spl_cluster) |>
    dplyr::count(spl_cluster) |>
    dplyr::arrange(spl_cluster) |>
    dplyr::pull(n)

  # Compute separator positions
  block_positions <- 0.5 + cumsum(block_sizes)[-length(block_sizes)]

  # Compute block sizes AFTER reordering (ASV)
  b_sizes <- df_long |>
    dplyr::distinct(taxon, tax_cluster) |>
    dplyr::count(tax_cluster) |>
    dplyr::arrange(tax_cluster) |>
    dplyr::pull(n)

  # Compute separator positions
  asv_positions <- 0.5 + cumsum(b_sizes)[-length(b_sizes)]

  mat <- ggplot(df_long, aes(
    x = Sample,
    y = taxon,
    fill = color_bin
  )) +
    coord_flip() +
    geom_tile() +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        face = "italic"
      ),
      axis.text.y = element_text()
    ) +
    scale_fill_manual(values = heatmap_color_bin) +
    scale_y_discrete(name = "ASV", limits = correct_taxon_order) +
    scale_x_discrete(name = "Sample", limits = correct_sample_order) +
    geom_vline(
      xintercept = block_positions,
      color = "black",
      linewidth = 2
    ) +
    geom_hline(
      yintercept = asv_positions,
      color = "black",
      linewidth = 2
    ) + labs(fill = legend_name)


  return(mat)
}
