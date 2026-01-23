#' @title Physeq filtering info
#'
#' @description gives the info of how many OTU and samples were lost before cleaning
#'
#' @param ps_before
#' @param ps_after
#'
#'
#' @return
#'
#' @export
sum_physeq_clean <- function(ps_before, ps_after) {
  reads_before <- sum(sample_sums(ps_before))
  reads_after  <- sum(sample_sums(ps_after))

  tibble::tibble(
    samples_before = nsamples(ps_before),
    samples_after  = nsamples(ps_after),
    samples_lost   = nsamples(ps_before) - nsamples(ps_after),
    samples_lost_pct =
      100 * (nsamples(ps_before) - nsamples(ps_after)) / nsamples(ps_before),

    OTUs_before    = ntaxa(ps_before),
    OTUs_after     = ntaxa(ps_after),
    OTUs_lost      = ntaxa(ps_before) - ntaxa(ps_after),
    OTUs_lost_pct  =
      100 * (ntaxa(ps_before) - ntaxa(ps_after)) / ntaxa(ps_before),

    reads_before   = reads_before,
    reads_after    = reads_after,
    reads_lost     = reads_before - reads_after,
    reads_lost_pct =
      100 * (reads_before - reads_after) / reads_before
  )
}
