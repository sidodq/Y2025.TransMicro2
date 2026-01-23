#' @title Add the marker at the end of the sample ID to distinguish 16S replicates from ITS
#'
#' @description
#'filter phyloseq object by suppressing ASVs with less than 70 % of prevalence per SampleType x Fruit
#'
#' @param ps phyloseq object
#' @param prev_thresh = 0.7
#' @param group_cols = c("SampleType", "FruitType")
#'
#'
#' @export
filter_phyloseq_prevalence <- function(physeq, group_cols = c("SampleType", "FruitType"), prevalence = 0.7) {

  stopifnot(inherits(physeq, "phyloseq"))

  # Sample data
  sdata <- as.data.frame(phyloseq::sample_data(physeq))
  sdata$Sample <- rownames(sdata)

  # OTU table
  otu <- as(phyloseq::otu_table(physeq), "matrix")
  if (!phyloseq::taxa_are_rows(physeq)) otu <- t(otu)

  otu_df <- as.data.frame(t(otu))
  otu_df$Sample <- rownames(otu_df)

  # Join metadata
  otu_long <- dplyr::left_join(
    otu_df,
    sdata[, c("Sample", group_cols)],
    by = "Sample"
  )

  # Presence / absence
  otu_long[, rownames(otu)] <- (otu_long[, rownames(otu)] > 0) * 1

  # Prevalence per group
  prev <- otu_long |>
    tidyr::pivot_longer(
      cols = rownames(otu),
      names_to = "OTU",
      values_to = "present"
    ) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), OTU) |>
    dplyr::summarise(prev = mean(present), .groups = "drop") |>
    dplyr::filter(prev >= prevalence)

  taxa_keep <- unique(prev$OTU)

  phyloseq::prune_taxa(taxa_keep, physeq)
}

#'
#'
# filter_phyloseq_by_prevalence <- function(ps, prev_thresh = 0.7, group_cols = c("SampleType", "FruitType")) {
#
#   sdata <- sample_data(ps) %>% as.data.frame()
#   otu <- as(otu_table(ps), "matrix")
#   if (!taxa_are_rows(ps)) otu <- t(otu)
#
#   otu_df <- as.data.frame(t(otu))
#   otu_df$Sample <- rownames(otu_df)
#   otu_df <- left_join(otu_df, sdata[, c("Sample", group_cols)], by = "Sample")
#
#   prev_tbl <- otu_df %>%
#     tidyr::pivot_longer(cols = -c(Sample, all_of(group_cols)),
#                         names_to = "OTU", values_to = "Abundance") %>%
#     mutate(present = Abundance > 0) %>%
#     group_by(across(all_of(group_cols)), OTU) %>%
#     summarise(prevalence = mean(present), .groups = "drop") %>%
#     filter(prevalence >= prev_thresh)
#
#   OTUs_to_keep <- unique(prev_tbl$OTU)
#
#   ps_filtered <- prune_taxa(OTUs_to_keep, ps)
#   ps_filtered <- prune_samples(sample_sums(ps_filtered) > 0, ps_filtered)
#   return(ps_filtered)
# }
