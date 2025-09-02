#' @title Create prevalence/relative abundance plot
#'
#' @param phyloseq_object the phyloseq object
#' @param taxonomic_rank "Family" or "Genus"
#' @param n_families max number of families to show
#' @param grouping_factors category on which the prevalence will be calculated ex. "LifeStage"
#' @param color_set set of colors to use
#'
#' @description
#' Plot the the prevalence of the species present -at chosen taxonomic level- in function of their relative abundance.
#'
#'
#' @export
Make_prevalence_abundance_plot_3 <- function(phyloseq_object, taxonomic_rank, n_families, grouping_factors, color_set) {

# Agglomerate to Family level (i.e. ASV from same family are summed per sample)
phyloseq_object <- phyloseq::tax_glom(phyloseq_object, taxrank = taxonomic_rank)

## Compute prevalence, relative abundance and specificity for each level of factors
prevalence_relativabundance <- phyloseq_object |>
  phyloseq.extended::estimate_prevalence(group = grouping_factors ) |>
  dplyr::filter(group != "na")

unique(prevalence_relativabundance$group)


## Add family for each cluster
tax_df <- as.data.frame(phyloseq::tax_table(phyloseq_object))
tax_df$otu <- rownames(tax_df) # OTUs as ID

prevalence_relativabundance <- dplyr::left_join(prevalence_relativabundance, tax_df, by = "otu")

# Select top N families by mean abundance
top_families <- prevalence_relativabundance |>
  dplyr::group_by(.data[[taxonomic_rank]]) |>
  dplyr::summarise(mean_abund = mean(abundance, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(desc(mean_abund)) |>
  dplyr::slice_head(n = n_families) |>
  dplyr::pull(.data[[taxonomic_rank]])
families_to_label <- c("Acetobacteraceae", "Erwiniaceae", "Lactobacillaceae","Enterobacteriaceae", "Bacillaceae", "Unknown", "Pichiaceae", "Metschnikowiaceae", "Saccharomycodaceae", "Saccharomycetaceae")

# Optional: manually label a few more if desired
families_to_label <- union(
  c("Acetobacteraceae", "Erwiniaceae", "Lactobacillaceae", "Enterobacteriaceae", "Bacillaceae", "Unknown",
    "Pichiaceae", "Metschnikowiaceae", "Saccharomycodaceae", "Saccharomycetaceae"),
  top_families)

p <- prevalence_relativabundance |>
  ggplot(aes(x = abundance, y = prevalence, color = .data[[taxonomic_rank]])) +
  geom_point(size = 3) +
  geom_text_repel(
    data = dplyr::filter(prevalence_relativabundance, .data[[taxonomic_rank]] %in% families_to_label),
    aes(label = .data[[taxonomic_rank]])
  ) +
  scale_color_manual(values = color_set) +
  theme_minimal() +
  labs(x = "Mean relative Abundance", y = "Prevalence") +
  facet_wrap(vars(group)) +
  theme(
    strip.text = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    axis.title.x = element_text(size = 13, face = "bold"),
    axis.title.y = element_text(size = 13, face = "bold")
  )

return(p)
}
