#' @title Subset a phyloseq object
#'
#' @description
#' subset a phyloseq object based on a condition to verify
#'
#' @param phyloseq_object the phyloseq object to subset
#' @param condition_expr condition under which the selection will be done ex. grepl("Focal_Sites", Sample_Category)
#'
#' @examples
#' S_ps_16S_rare_FocalSites_All <- subset_samples(S_ps_rare_16S, grepl("Focal_Sites", Sample_Category)) # not only Sample_Category == "All_stages_focal_Sites"
#' S_ps_16S_rare_FocalSites_All <- prune_taxa(taxa_sums(S_ps_16S_rare_FocalSites_All)>0, S_ps_16S_rare_FocalSites_All)
#'
#' @return
#' the subseted phyloseq object
#'
#' @export
Subset_physeq <- function(phyloseq_object, condition_expr) {
  condition_expr <- rlang::enquo(condition_expr)

  # Convert sample data to data frame for evaluation
  sample_df <- as.data.frame(sample_data(phyloseq_object))

  # Evaluate condition
  keep_samples <- rlang::eval_tidy(condition_expr, data = sample_df)

  # Subset phyloseq object
  subset_physeq_obj <- prune_samples(keep_samples, phyloseq_object)
  subset_physeq_obj <- prune_taxa(taxa_sums(subset_physeq_obj) > 0, subset_physeq_obj)

  return(subset_physeq_obj)
}
