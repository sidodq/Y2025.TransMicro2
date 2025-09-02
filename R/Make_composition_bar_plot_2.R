#' @title Make composition bar plot 2
#'
#' @description
#' Function to make bar plot for bacterial or fungal taxa at the family/genus level, regrouping samples by groups
#'
#'
#' @param phyloseq_object_filtered the phyloseq object
#' @param kingdom the kingdom of interest : "Fungi" or "Bacteria"
#' @param taxonomic_level the variable to be represented by the colours in the graph, here the taxa level : "Family" or "Genus"
#' @param TaxaNb the number of taxa to show (min 10)
#' @param nested_modality the modality(ies) under which the groups will be made ex. ~LifeStage
#' @param color_set the set of colors to use ex. color_16S
#'
#'
#' @return
#' The bar plot with in x the group of samples tested and in y the abundance of taxa families/genus.
#'
#' @export
Make_composition_bar_plot_2 <- function(phyloseq_object_filtered, kingdom, taxonomic_level, TaxaNb, nested_modality, color_set) {

plot <- phyloseq.extended::plot_composition(
  phyloseq_object_filtered,
  "Kingdom", kingdom,
  taxonomic_level,
  numberOfTaxa = TaxaNb,
  fill = taxonomic_level
) +
  ggh4x::facet_nested(nested_modality, scales = "free_x", space = "fixed") +
  ggplot2::scale_fill_manual(values = color_set) +
  ggplot2::labs(x = "Samples", y = "Relative abundance", title = "") +
  ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
  ggplot2::theme(
    strip.text.x = ggplot2::element_text(size = 13, face = "bold"),
    axis.text.x = ggplot2::element_text(size = 8),
    axis.text.y = ggplot2::element_text(size = 12),
    axis.title.x = ggplot2::element_text(size = 13),
    axis.title.y = ggplot2::element_text(size = 13),
    legend.text = ggplot2::element_text(size = 11),
    legend.title = ggplot2::element_text(size = 13))
plot(plot)
}
