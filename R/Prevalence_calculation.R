#' @title Calculate the prevalence
#'
#' @description
#' calculate the prevalence of clusters within a given group of samples
#'
#' @param physeq_object physeq object
#' @param metadata contains all the variables susceptible to be used for calculating the prevalence on
#' @param group the reference variable (group of sample) for which the prevalence will be calculated
#
#'
#' @return a table with the prevalence calculated per cluster (+the specificity and the abundance) compared to the bigger group
#'
#'
#' @export
prevalence_by <- function(physeq_object, metadata, group) {
  group_sym <- rlang::ensym(group)  # handles both quoted and unquoted

  metadata |>
    dplyr::group_by(!!group_sym) |>
    dplyr::summarise(prevalence = estimate_prevalence(physeq_object, group = !!group_sym))
}

