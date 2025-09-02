#' @title Make composition bar plot
#'
#' @description
#' Function to make bar plot for bacterial or fungal taxa at the family/genus level.
#'
#'
#' @param phyloseq_object the phyloseq object
#' @param x_variable the variable included in the phyloseq object, which will correspond to the x axis of the graph
#' @param taxonomic_level the variable to be represented by the colours in the graph, here the taxa level : "Family" or "Genus"
#' @param abundance_modality ="fill" for normalized abundance, = "stack" for absolute abundance
#'
#'
#' @return
#' The bar plot with in x the group of samples tested and in y the abundance of taxa families/genus.
#'
#' @export
Make_composition_bar_plot <- function(phyloseq_object, x_variable, taxonomic_level, abundance_modality) {
  ggplot(phyloseq_object, aes_string(x = x_variable, y = "Abundance", fill = taxonomic_level)) +
    geom_bar(stat = "identity", color = "black", position = abundance_modality)
}
