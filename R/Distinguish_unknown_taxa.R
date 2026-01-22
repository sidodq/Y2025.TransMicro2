#' @title Distinguishing unknown taxa
#'
#' @description
#' Number the "Unknown" taxa corresponding with clusters' numbers
#'
#' @param phyloseq_object the phyloseq object to subset
#' @param ranks taxo ranks to treat like "ranks = c("Family", "Genus")"
#' @param outfile file direction where to save the new taxonomic table
#'
#'
#' @export
fix_unknown_taxa <- function(phyloseq_object, ranks, outfile) {
  tax <- as.data.frame(tax_table(phyloseq_object))
  otu_names <- rownames(tax)

  for (rank in ranks) {
    tax[[rank]] <- ifelse(tax[[rank]] == "Unknown", paste0("Unknown_", sub("Cluster_", "", otu_names)), as.character(tax[[rank]]))
  }

  tax_table(phyloseq_object) <- as.matrix(tax)

  write.table(tax_table(phyloseq_object), file = outfile, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

  return(phyloseq_object)
}

# tax <- as.data.frame(tax_table(ps_rare_16S)) # Extract taxonomy table as a data frame
# otu_names <- rownames(tax) # Extract OTU names
#
# tax$Family <- ifelse(
#   tax$Family == "Unknown",
#   paste0("Unknown_", sub("Cluster_", "", otu_names)), # Add suffix with cluster number if Family == "Unknown"
#   as.character(tax$Family))
#
#
#
# tax_table(ps_rare_16S) <- as.matrix(tax) # Put it back into the phyloseq object
# write.table(tax_table(ps_rare_16S), file=file.path(here::here(), "data", "derived_data", "3_CommunityComposition", "16S_TM123_tax_table_rare.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#
#
# tax <- as.data.frame(tax_table(ps_rare_ITS)) # Extract taxonomy table as a data frame
# otu_names <- rownames(tax) # Extract OTU names
#
# tax$Family <- ifelse(
#   tax$Family == "Unknown",
#   paste0("Unknown_", sub("Cluster_", "", otu_names)), # Add suffix with cluster number if Family == "Unknown"
#   as.character(tax$Family))
#
# tax_table(ps_rare_ITS) <- as.matrix(tax) # Put it back into the phyloseq object
# write.table(tax_table(ps_rare_ITS), file=file.path(here::here(), "data", "derived_data", "3_CommunityComposition", "ITS_TM123_tax_table_rare.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
# ```
#
# ``` {r}
# # Distinguish "Unknown" GENUS by numeroting them with their corresponding cluster number
#
# tax <- as.data.frame(tax_table(ps_rare_16S)) # Extract taxonomy table as a data frame
# otu_names <- rownames(tax) # Extract OTU names
#
# tax$Genus <- ifelse(
#   tax$Genus == "Unknown",
#   paste0("Unknown_", sub("Cluster_", "", otu_names)), # Add suffix with cluster number if Family == "Unknown"
#   as.character(tax$Genus))
#
# tax_table(ps_rare_16S) <- as.matrix(tax) # Put it back into the phyloseq object
# write.table(tax_table(ps_rare_16S), file=file.path(here::here(), "data", "derived_data", "3_CommunityComposition", "16S_TM123_tax_table_rare.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#
#
# tax <- as.data.frame(tax_table(ps_rare_ITS)) # Extract taxonomy table as a data frame
# otu_names <- rownames(tax) # Extract OTU names
#
# tax$Genus <- ifelse(
#   tax$Genus == "Unknown",
#   paste0("Unknown_", sub("Cluster_", "", otu_names)), # Add suffix with cluster number if Family == "Unknown"
#   as.character(tax$Genus))
#
# tax_table(ps_rare_ITS) <- as.matrix(tax) # Put it back into the phyloseq object
#write.table(tax_table(ps_rare_ITS), file=file.path(here::here(), "data", "derived_data", "3_CommunityComposition", "ITS_TM123_tax_table_rare.tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
