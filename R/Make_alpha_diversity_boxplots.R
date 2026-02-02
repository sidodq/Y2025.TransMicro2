#' @title Create alpha diversity boxplots
#'
#' @param phyloseq_object the phyloseq object
#' @param index the alpha diversity index
#' @param x_variable the variable included in the phyloseq object, which will correspond to the x axis of the graph
#' @param fill the variable to be represented by the colours in the graph
#' @param fill_color a character string with specific colours
#' @param name_fill_modality a character string with specific names of the fill variable
#' @param variable_separated1 a variable for separating the boxplots
#' @param variable_separated2 another variable for separating the boxplots
#' @param y_lab name of the y-axis
#'
#' @description
#' Make a boxplot for alpha diversity index, and define the desired explicative variables.
#'
#' @return
#' A figure for alpha diversity index. It's a "ggplot" class object.
#'
#' @examples
#' divalpha_ML_MS_G0W_G0GF_G1_G2_Pop_16S <- Make_alpha_diversity_boxplots(metadata_alphadiv = metadata_ML_MS_G0W_G0GF_G1_G2_16S_alphadiv, x_axis = "Site_Short", fill = "Site_Short", color = Population)
#'
#' @export
Make_alpha_diversity_boxplots <- function(phyloseq_object, index, x_variable, y_lab, fill, fill_color = NULL, name_fill_modality = NULL, variable_separated1 = NULL, variable_separated2 = NULL)
  {{

  alpharich <- estimate_richness(phyloseq_object, measures=c(index))
  metadata_alphadiv <- data.frame(alpharich, sample_data(phyloseq_object), check.names=TRUE)

  Figure <- ggplot(metadata_alphadiv, aes(x = metadata_alphadiv[,x_variable], y = metadata_alphadiv[,index] , fill = metadata_alphadiv[,fill]))+
    geom_boxplot(outlier.shape = NA) +
    labs(x="", y = y_lab, fill = fill) +
    theme(axis.text.x = element_text(face=, color="Black", size=18, angle=90), axis.text.y = element_text(face=, color="Black", size=18, angle=0), legend.position = "right") +
    theme(axis.title.x = element_text(color="Black", size=18, face="bold"),
          axis.title.y = element_text(color="Black", size=18, face="bold"),
          legend.title = element_text(color = "Black", size = 18, face="bold"),
          legend.text = element_text(color = "Black", size = 18)) +
    scale_fill_manual(values =  fill_color, labels = name_fill_modality)+
    geom_jitter(size = 1.5, shape = 21, position = position_jitter(0.1), alpha = 0.9, color = "black")+
    scale_x_discrete(labels = NULL)

  if (!is.null(variable_separated1) & !is.null(variable_separated2)) {
    Figure <- Figure + ggh4x::facet_nested(~as.character(metadata_alphadiv[,variable_separated1]) + as.character(metadata_alphadiv[,variable_separated2]), scales = "free_x", space = "free_x") +
      theme_bw()}


  }
return(list(Figure = Figure, metadata_alphadiv = metadata_alphadiv))

}
