#' @title Create prevalence/abundance plot
#'
#' @param plot_df the data frame with the data to plot (=prevalence calculated by ASVs)
#' @param taxonomic_rank "Family" or "Genus"
#' @param color_set the set of colors to use
#' @param title_plot clear enough
#' @param relative TRUE if relative abundance, otherwise FALSE
#'
#' @description
#' Plot the the prevalence of the species present -at chosen taxonomic level- in function of their abundance
#'
#' @return
#' NA
#'
#' @export
Make_prevalence_abundance_plot2 <- function(plot_df, taxonomic_rank, color_set, title_plot, relative) {

  plot_df <- plot_df %>%
    group_by(GroupingVariable, group, .data[[taxonomic_rank]]) %>%  # group: your sample group, Family from taxonomy
    dplyr::summarise(
      prevalence = mean(prevalence, na.rm = TRUE),
      abundance = sum(abundance, na.rm = TRUE),
      specificity = mean(specificity, na.rm = TRUE), # if you want to keep it
      .groups = "drop")

  p <- ggplot(plot_df, aes(x = abundance, y = prevalence, label = .data[[taxonomic_rank]])) +
      geom_point(aes(color = .data[[taxonomic_rank]]), size = 3) +
      geom_text_repel(max.overlaps = 20) +
      scale_color_manual(values = color_set) +
      theme_minimal() +
      labs(x = "Abundance", y = "Prevalence", title = title_plot)

    if (relative) {
      p <- p + facet_wrap(~GroupingVariable)}

    return(p)
}
