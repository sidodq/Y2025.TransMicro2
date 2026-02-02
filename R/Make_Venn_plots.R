#' @title Extract the taxa, list of Families or Genera, per tested variable (FruitType, LifeStage, Pop, etc.)
#'
#' @description Function to extract genera for a given category
#'
#' @param physeq the phyloseq object
#' @param variable the sample category of interest (ex. FruitType), to write between ""
#' @param value the value of the category (ex. BL), to write between ""
#'
#'
#' @export
get_otus_for_type <- function(physeq, variable, value) {
# Extract metadata
meta <- as.data.frame(sample_data(physeq))

# Subset manually
sample_keep <- rownames(meta)[meta[[variable]] == value]
ps_sub <- prune_samples(sample_keep, physeq)
ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)

# Return ASV names
return(taxa_names(ps_sub))
}

# @param taxo the taxonomic level wanted (ex. Family)
# get_taxa_for_type <- function(physeq, variable, value, taxo) {
#   # Subset manually using standard R logic
#   sample_keep <- sample_names(physeq)[get_variable(physeq, variable) == value]
#   ps_sub <- prune_samples(sample_keep, physeq)
#   ps_sub <- prune_taxa(taxa_sums(ps_sub) > 0, ps_sub)
#   tax_df <- as.data.frame(tax_table(ps_sub))
#   unique(tax_df[[taxo]])
# }


#' @title Make Venn plot divided per FruitType (BL, CH, ST)
#'
#' @description
#' Function to represent Families or Genus shared between and accross Fruit Types
#'
#'
#' @param list_BL list of the Families or Genus proper to BL samples selected
#' @param list_CH idem with samples from CH
#' @param list_ST idem with samples from ST
#'
#'
#' @export
Make_Venn_plot_FruitType <- function(list_BL, list_CH, list_ST) {
  # Basic 3-set Venn Diagram
  venn.plot <- draw.triple.venn(
    area1 = length(list_BL), area2 = length(list_CH), area3 = length(list_ST),
    n12 = length(intersect(list_BL, list_CH)), n23 = length(intersect(list_CH, list_ST)), n13 = length(intersect(list_BL, list_ST)),
    n123 = length(Reduce(intersect, list(list_BL, list_CH, list_ST))),
    category = c("Blackberry (BL)", "Cherry (CH)", "Strawberry (ST)"),
    fill = c("purple4", "red4", "red2"),
    cat.cex = 1.2, cex = 1.3, fontfamily = "sans", cat.pos = c(-20, 20, 0), cat.dist = c(0.06, 0.06, - 0.45), scaled = TRUE)

  # Draw the diagram
  grid::grid.newpage()  # Start fresh for the plot
  grid::grid.draw(venn.plot)

  return(venn.plot)
}

#' @title Make Venn plot divided per LifeStage (Mother, Larvae, Pupae, PupalCase, Daughter)
#'
#' @description
#' Function to represent Families or Genus shared between and across Life Stages
#'
#'
#' @param list_THU list of the Families or Genus proper to Mothers samples
#' @param list_FIN idem with samples from Larvae
#' @param list_LOM idem with samples from Pupae
#' @param list_IST idem with samples from Pupal Cases
#' @param list_MIL idem with samples from Daughters
#' @param list_VAR
#'
#'
#' @export
Make_Venn_plot_FruitSite <- function(list_THU, list_FIN, list_LOM, list_IST, list_MIL, list_VAR) {

  venn_data <- list(THURINS = list_THU, FINISTRET = list_FIN, LOSMASOS = list_LOM, ILESURTET = list_IST, MILLAS = list_MIL, VARETZ = list_VAR)

ggVennDiagram(venn_data, label_alpha = 0, edge_size = 0.8) +
  scale_fill_gradient(low = "white", high = "#8B0A50", name ="ASVs number") + theme_void() + theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10))
}



#' @title Make Venn plot divided per LifeStage (Mother, Larvae, Pupae, PupalCase, Daughter)
#'
#' @description
#' Function to represent Families or Genus shared between and across Life Stages
#'
#'
#' @param list_M list of the Families or Genus proper to Mothers samples
#' @param list_L idem with samples from Larvae
#' @param list_P idem with samples from Pupae
#' @param list_PC idem with samples from Pupal Cases
#' @param list_D idem with samples from Daughters
#'
#'
#' @export
Make_Venn_plot_LifeStage <- function(list_M, list_L, list_P, list_PC, list_D) {

  venn_data <- list(Mother = list_M, Larva = list_L, Pupa = list_P, PupalCase = list_PC, Daughter = list_D)

  ggVennDiagram(venn_data, label_alpha = 0, edge_size = 0.8) +
    scale_fill_gradient(low = "white", high = "steelblue", name ="ASVs number") + theme_void() + theme(legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 10))
}





#' @title Common ASVs
#'
#' @description extract the list of ASVs shared by a % of samples
#'
#'
#' @param physeq phyloseq object
#' @param thresh percentage of samples
#'
#'
#' @export
Common_ASVs <- function(physeq, thresh) {

# Extract the abundance matrix
otu_mat <- as(otu_table(physeq), "matrix")

# If taxa are in rows, nothing to change; if in columns, transpose it
if (taxa_are_rows(physeq)) {
  otu_mat <- otu_mat
} else {
  otu_mat <- t(otu_mat)}

# Compute presence/absence matrix
otu_pa <- otu_mat > 0

# Find taxa present in at least 80% of samples
threshold <- thresh * ncol(otu_pa)
shared_ASVs <- rownames(otu_pa)[rowSums(otu_pa) >= threshold]

return(shared_ASVs)
}



#'
# Make_Venn_plot_LifeStage <- function(list_M, list_L, list_P, list_PC, list_D) {
#
#   venn.plot <- draw.quintuple.venn(
#     area1 = length(list_M),
#     area2 = length(list_L),
#     area3 = length(list_P),
#     area4 = length(list_PC),
#     area5 = length(list_D),
#
#     n12 = length(intersect(list_M, list_L)),
#     n13 = length(intersect(list_M, list_P)),
#     n14 = length(intersect(list_M, list_PC)),
#     n15 = length(intersect(list_M, list_D)),
#     n23 = length(intersect(list_L, list_P)),
#     n24 = length(intersect(list_L, list_PC)),
#     n25 = length(intersect(list_L, list_D)),
#     n34 = length(intersect(list_P, list_PC)),
#     n35 = length(intersect(list_P, list_D)),
#     n45 = length(intersect(list_PC, list_D)),
#
#     # 3-way intersections
#     n123 = length(Reduce(intersect, list(list_M, list_L, list_P))),
#     n124 = length(Reduce(intersect, list(list_M, list_L, list_PC))),
#     n125 = length(Reduce(intersect, list(list_M, list_L, list_D))),
#     n134 = length(Reduce(intersect, list(list_M, list_P, list_PC))),
#     n135 = length(Reduce(intersect, list(list_M, list_P, list_D))),
#     n145 = length(Reduce(intersect, list(list_M, list_PC, list_D))),
#     n234 = length(Reduce(intersect, list(list_L, list_P, list_PC))),
#     n235 = length(Reduce(intersect, list(list_L, list_P, list_D))),
#     n245 = length(Reduce(intersect, list(list_L, list_PC, list_D))),
#     n345 = length(Reduce(intersect, list(list_P, list_PC, list_D))),
#
#     # 4-way intersections
#     n1234 = length(Reduce(intersect, list(list_M, list_L, list_P, list_PC))),
#     n1235 = length(Reduce(intersect, list(list_M, list_L, list_P, list_D))),
#     n1245 = length(Reduce(intersect, list(list_M, list_L, list_PC, list_D))),
#     n1345 = length(Reduce(intersect, list(list_M, list_P, list_PC, list_D))),
#     n2345 = length(Reduce(intersect, list(list_L, list_P, list_PC, list_D))),
#
#     # 5-way intersection
#     n12345 = length(Reduce(intersect, list(list_M, list_L, list_P, list_PC, list_D))),
#
#     category = c("Mother (M)", "Larva (L)", "Pupa (P)", "Pupal Case (PC)", "Daughter (D)"),
#     fill = c("darkolivegreen3", "green3", "darkorange2", "#8B4500", "goldenrod2"),
#     cat.cex = 1.1,
#     cex = 1.2,
#     fontfamily = "sans",
#     scaled = TRUE
#   )
#
#   grid::grid.newpage()
#   grid::grid.draw(venn.plot)
#
#   return(venn.plot)
# }
