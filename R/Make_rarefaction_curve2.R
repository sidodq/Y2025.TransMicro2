#' @title Plot rarefaction curves
#' @param phylo_object a phyloseq object
#' @param x_intersection_line a number for positioning the dotted vertical line
#' @param color_variable a variable that distinguishes its modalities by different colors
#' @param maxdepth maximum depth for the x-axis of the graph
#'
#' @description
#' Create a rarefaction curve from a phyloseq object.
#'
#'
#' @export
Make_rarefaction_curve2 <- function(phylo_object, x_intersection_line, color_variable, title, maxdepth = FALSE, logy = FALSE) {

rarefaction_curve <- phyloseq.extended::ggrare(
  physeq = phylo_object,
  step = 50,
  color = color_variable,
  se = FALSE
) +
  theme_bw() +
  geom_vline(xintercept = x_intersection_line, color = "black", linetype = "dashed") +
  ggtitle(title)

if (maxdepth > 0) {
  rarefaction_curve <- rarefaction_curve + xlim(0, maxdepth)
}

if (logy) {
  rarefaction_curve <- rarefaction_curve + scale_y_log10()
}

return(rarefaction_curve)
}

# Make_rarefaction_curve2 <- function(phylo_object, x_intersection_line, color_variable, title, maxdepth=FALSE) {
#
#   rarefaction_curve <- phyloseq.extended::ggrare(physeq = phylo_object, step = 50, color = color_variable, se = FALSE) +
#     theme_bw() +
#     geom_vline(xintercept = x_intersection_line, color = "black", linetype = "dashed")+
#     ggtitle(title)
#
#   if(maxdepth>0){
#     rarefaction_curve <- rarefaction_curve + xlim(0, maxdepth)
#   }
#   return(rarefaction_curve)
#
# }
