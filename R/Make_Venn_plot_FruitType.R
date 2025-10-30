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
