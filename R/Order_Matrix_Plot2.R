#' @title Create heatmap with LBM clusters
#'
#' @param SBMobject BipartiteSBM object
#' @param BF_mat binary matrix
#' @param variable the variable included in the SBM object
#' @param varibale_name name of the variable for colouring the squares
#' @param modality_string modalities of the variable to colour the squares displayed in the legend
#' @param color_string character string with colours for each mode of the variable
#'
#' @description
#' Creation of a heatmap with samples in y and OTUs in x. A white square represents the absence of the OTU in the sample. A coloured square represents the presence of the OTU in the sample. Samples and OTUs are grouped according to the membership provided by the Latent Block Model. The colour of the squares is defined according to an environmental variable that can be chosen.
#'
#' @return
#' It's a "ggplot" class object.
#'
#' @examples
#' LBM_16S_MPC3_Fruit <- Order_Matrix_Plot(SBMobject = LBM_bernoulli_MPC3_16S, BF_mat = BF_MPC3_16S, color_variable = "Fruit", varibale_name = "Fruit", modality_string = c("Blackberry", "Cherry", "Strawberry"), color_string = Color_Fruit)

#'
#' @export
Order_Matrix_Plot2 <- function(SBMobject, BF_mat, mtd, color_variable, modality_string, variable_name, color_string, new_rows) {

  df = as.data.frame(BF_mat)
  df$taxon = rownames(df)
  df_long <-
    tidyr::pivot_longer(df,
                        cols = names(df)[-(dim(BF_mat)[[2]]+1)],
                        names_to = "Sample",
                        values_to = "PresenceAbsence")
  df_long <- dplyr::left_join(df_long, mtd, by = "Sample")

  Qrow <- SBMobject$nbBlocks["row"][[1]]
  Qcol <- SBMobject$nbBlocks["col"][[1]]

  df_tax_cluster = data.frame(taxon = rownames(BF_mat),
                              tax_cluster = SBMobject$memberships$col)
  df_long <- dplyr::left_join(df_long, df_tax_cluster, by = "taxon")

  df_long <- df_long %>% dplyr::arrange(tax_cluster)
  correct_taxon_order <- unique(df_long$taxon)

  df_spl_cluster = data.frame(Sample = colnames(BF_mat),
                              spl_cluster = new_rows)
  df_long <- dplyr::left_join(df_long, df_spl_cluster, by = "Sample")

  df_long <- df_long %>% dplyr::arrange(spl_cluster)
  correct_sample_order <- unique(df_long$Sample)
  print(df_long$FruitType[df_long$FruitType!="BL"])
  df_long[[color_variable]] <- factor(df_long[[color_variable]])

  df_long$color_bin <- factor(ifelse(df_long$PresenceAbsence==0, 0, paste(df_long$PresenceAbsence, df_long[[color_variable]], sep = "_")))
  levels(df_long$color_bin) <- c("", modality_string)
  heatmap_color_bin <- c("white", color_string)

  df_long$PresenceAbsence <- as.factor(df_long$PresenceAbsence)

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
      xintercept = 0.5 + cumsum(table(new_rows))[-length(table(new_rows))],
      color = "black",
      linewidth = 2
    ) +
    geom_hline(
      yintercept = 0.5 + cumsum(table(SBMobject$memberships$col))[-length(table(SBMobject$memberships$col))],
      color = "black",
      linewidth = 2
    ) + labs(fill = variable_name)


  return(mat)
}
