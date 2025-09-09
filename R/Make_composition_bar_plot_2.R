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
#' @export
Make_composition_bar_plot_2 <- function(phyloseq_object_filtered, kingdom, taxonomic_level, TaxaNb, nested_modality, color_set) {

  # Melt phyloseq object into tidy data
  df <- phyloseq::psmelt(phyloseq_object_filtered)

  # Filter to kingdom of interest
  df <- df[df$Kingdom == kingdom, ]

  # Collapse all "Unknown*" into a single "Unknown"
  df[[taxonomic_level]] <- ifelse(grepl("^Unknown(_[0-9]+)?$", df[[taxonomic_level]]),
                                  "Unknown",
                                  as.character(df[[taxonomic_level]]))

  # Keep only the top N taxa, the rest → "Other"
  top_taxa <- names(sort(tapply(df$Abundance, df[[taxonomic_level]], sum), decreasing = TRUE))[1:TaxaNb]
  df[[taxonomic_level]] <- ifelse(df[[taxonomic_level]] %in% top_taxa, df[[taxonomic_level]], "Other")

  # Force legend order: all taxa, then Unknown, then Other
  taxa_levels <- c(setdiff(unique(df[[taxonomic_level]]), c("Unknown", "Other")), "Unknown", "Other")
  df[[taxonomic_level]] <- factor(df[[taxonomic_level]], levels = taxa_levels)

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
