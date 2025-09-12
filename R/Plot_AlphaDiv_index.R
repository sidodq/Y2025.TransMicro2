#' @title plot alpha diversity index(es)
#'
#' @param df a dataframe obtained from the phyloseq object studied
#' @param group
#' @param group_labels
#' @param index
#' @param color_set
#'
#'
#' @export
plot_AlphaDiv_index <- function(df, group, group_labels, index, color_set) {
df <- df |>
  dplyr::filter((.data[[group]]) !="na") |>
  dplyr::mutate(
    !!group := factor(.data[[group]], levels = group_labels))

ggplot(df, aes(x = .data[[group]], y = .data[[index]], fill = .data[[group]])) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(1)) +
  geom_jitter(size = 1.5, shape = 21, position = position_jitter(0.1), alpha = 0.9, color = "black") +
  scale_fill_manual(values = color_set, labels = group_labels) +
  labs(x = group, y = paste(index)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "Black", size = 14, angle = 90),
    axis.text.y = element_text(color = "Black", size = 14),
    axis.title.x = element_text(color = "Black", size = 14, face = "bold"),
    axis.title.y = element_text(color = "Black", size = 14, face = "bold"),
    legend.title = element_text(color = "Black", size = 14, face = "bold"),
    legend.text = element_text(color = "Black", size = 14),
    legend.position = "right")
# +
#   # 👇 Add significance stars + required a new parameter "method"
#   stat_compare_means(comparisons = combn(group_labels, 2, simplify = FALSE), method = method,  label = "p.signif", tip.length = 0.01, size = 5)
}
