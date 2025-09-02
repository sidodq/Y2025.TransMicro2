#' @title Create prevalence/relative abundance plot
#'
#' @param phyloseq_object the phyloseq object
#' @param taxonomic_rank "Family" or "Genus"
#' @param grouping_factors category on which the prevalence will be calculated ex. "LifeStage"
#' @param color_set set of colors to use
#'
#' @description
#' Plot the the prevalence of the species present -at chosen taxonomic level- in function of their relative abundance.
#'
#'
#' @export
Make_prevalence_abundance_plot <- function(data, color_set) {


  families_to_label <- c("Acetobacteraceae", "Erwiniaceae", "Lactobacillaceae","Enterobacteriaceae", "Bacillaceae", "Unknown", "Pichiaceae", "Metschnikowiaceae", "Saccharomycodaceae", "Saccharomycetaceae")

  p <- data |>
    ggplot(aes(x = RelativeAbundance, y = Prevalence, color = Family)) +
    geom_point(size = 3) +

    geom_text_repel(
      data = dplyr::filter(data, Family %in% families_to_label),
      aes(label = Family)) +
    scale_color_manual(values = color_set) +
    theme_minimal() +
    labs(x = "Relative Abundance", y = "Prevalence") +
    theme(
      strip.text = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11), legend.title = element_text(size = 11, face = "bold"),
      legend.justification = c("left", "bottom"),
      axis.title.x = element_text(size = 13, face = "bold"), axis.title.y = element_text(size = 13, face = "bold"))

  return(p)
}
