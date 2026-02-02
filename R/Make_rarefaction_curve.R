#' @title Make a rarefaction curve
#' @param phylo_object a phyloseq object
#' @param x_intersection_line a number for positioning the dotted vertical line
#' @param variable a variable that can be graphed by modality of this variable
#' @param color_variable a variable that distinguishes its modalities by different colors
#' @param number_columns number to set the number of columns for graphics
#' @param maxdepth maximum depth for the x-axis of the graph
#'
#' @description
#' Create a rarefaction curve from a phyloseq object.
#'
#' @examples
#' rarefaction_curve_ML_MS_G0W_G0GF_G1_G2_16S <- Make_rarefaction_curve(phylo_object = PHYLO_ML_MS_G0W_G0GF_G1_G2_16S, x_intersection_line = 6000, variable = "Generation_SplitG0", number_columns = 4 , color_variable = "Generation_SplitG0")
#'
#' @return
#' A rarefaction curve
#'
#' @export
Make_rarefaction_curve <- function(phylo_object, x_intersection_line, variable, color_variable, number_columns, title, maxdepth=FALSE) {

  rarefaction_curve <- phyloseq.extended::ggrare(physeq = phylo_object, step = 50, color = color_variable, se = FALSE) +
    facet_wrap(as.formula(paste("~", variable)), ncol = number_columns) +
    theme_bw() +
    geom_vline(xintercept = x_intersection_line, color = "black", linetype = "dashed")+
    ggtitle(title)

  if(maxdepth>0){
    rarefaction_curve <- rarefaction_curve + xlim(0, maxdepth)
  }
  return(rarefaction_curve)

}

