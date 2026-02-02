#' @title plot alpha diversity index(es) on multiple sample levels (BoxPlot 1 index & 2 factors)
#'
#' @param df a dataframe obtained from the phyloseq object studied
#' @param group
#' @param group_labels
#' @param index
#' @param color_set
#' @param facet_factor description
#'
#'
#' @export
plot_AlphaDiv_index_groups <- function(df, group, group_labels, index, color_set, facet_factor) {
  df <- df |>
    dplyr::filter((.data[[group]]) !="na") |>
    dplyr::filter((.data[[facet_factor]]) !="na") |>
    dplyr::mutate(!!group := factor(.data[[group]], levels = group_labels))

  ggplot(df, aes(x = .data[[group]], y = .data[[index]] , fill = .data[[facet_factor]]))+
    geom_boxplot(outlier.shape = NA, position=position_dodge(1)) +
    theme_classic() +
    labs(x=group, y = paste(index)) +
    theme(axis.text.x = element_blank(), axis.text.y = element_text(face=, color="Black", size=14, angle=0), legend.position = "right") +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(color="Black", size=14, face="bold"),
          legend.title = element_text(color = "Black", size = 14, face="bold"),
          legend.text = element_text(color = "Black", size = 14)) +
    scale_fill_manual(values =  color_set, labels = c("BL","CH","ST"))+
    geom_jitter(size = 1.5, shape = 21, position = position_jitter(0.1), alpha = 0.9, color = "black") + ggh4x::facet_nested(as.formula(paste("~", group, "+", facet_factor)), scales = "free_x", space = "free_x")
}
