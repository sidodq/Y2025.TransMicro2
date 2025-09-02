#' @title Filter contaminant
#' @param x abundance table (only samples that need to be filtered)
#' @param Tcc_Tfa threshold (a number) to be subtracted from the abundances of each cluster for each sample
#' @description
#' Function that filters the abundance table :
#'- if the number of reads (=sequences) is less than the threshold then 0 is indicated
#'- otherwise the threshold is subtracted from the number of reads
#'
#' Threshold TCC data are based on the maximum number of sequences observed in a negative or positive control for a particular OTU in each run.
#' Threshold TFA data are based on the false-assignment rate (0.02%) weighted by the total number of sequences of each OTU in each run.
#' Galan et al. (2016)
#'
#' @examples
#'Rfa <- 0.0001
#'
#'at_MS_16S_ML_TccTfafiltered <- at_MS_16S_ML %>%
#'  dplyr::mutate(dplyr::across(
#'    -contains("CT-EXT") & -contains("CT-PCR") & -contains("NC-Index") &
#'      (ends_with("-A")|ends_with("-B")|ends_with("-C")|contains("INTERN")|contains("MOCK")|contains("Undetermined")),
#'   ~Filter_contaminant(., at_MS_16S_ML_Tccfiltered$observation_sum*Rfa)))
#'
#'
#' @return
#' An abundance table filtered by a threshold.
#'
#'
#' @export
#'
Filter_contaminant <- function(x, Tcc_Tfa){
  x <- ifelse(x <= Tcc_Tfa, 0, (x - Tcc_Tfa))
  return(x)
}

