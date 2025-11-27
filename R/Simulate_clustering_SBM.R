#' @title Compute clustering (memberships) on each simulated matrix with SBM algorithm
#' @param  NPERM number of permutations performed on the matrix
#' @param FB matrix with samples in rows and taxa in columns
#' @param name name of exported file
#' @description
#' This function is used to Compute clustering (memberships) on each simulated matrix with SBM algorithm.
#'
#' @return
#' A table with the memberships of each sample (in rows) as a function of each permuted matrix (in columns).
#' The number of columns corresponds to the number of permutations and the number of rows corresponds to the number of samples.
#'
#' @examples
#' simulated_clustering_MPC_16S_GENUS <- Simulate_Clustering_SBM(FB = FB_bin_MPC_GENUS_16S, NPERM = 1000, name = "simulated_clustering_MPC_16S_GENUS.tsv")
#'
#'
#' #' @export
#' Simulate_Clustering_SBM <- function(FB, NPERM, name)
#'
#' {
#'   # Produce NPERM simulated matrix
#'   .Random.seed
#'   netBin.config <- simulate(nullmodel(as.data.frame(FB), "curveball"), nsim = NPERM)
#'
#'   # Compute clustering on each simulated graph with SBM algorithm
#'   n_sim = 0
#'   list_sim = {}
#'
#'   while (n_sim < NPERM)
#'   {
#'     one_batch = sapply(1:NPERM, function(x) {
#'       simulSBM = estimateBipartiteSBM(as.matrix(netBin.config[, , x]),
#'                                       model = 'bernoulli',
#'                                       estimOptions = list(verbosity = 1, plot = F))
#'       c(simulSBM$memberships$row)
#'     })
#'     n_sim <- n_sim + NPERM
#'     list_sim <- cbind(list_sim, one_batch)
#'   }
#'
#'   simulated_clustering <- list_sim %>% as.data.frame()
#'
#'   # Export
#'   write.table(simulated_clustering, file=file.path(here::here(), "data", "derived_data", "5_SBMclustering", name), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#'
#'   return(simulated_clustering)
#' }




#' @title Compute simulated matrix with parallel simulations
#'
#' @param  NPERM number of permutations performed on the matrix
#' @param FB matrix with samples in rows and taxa in columns
#' @param name name of exported file
#' @description
#' This function is used to Compute clustering (memberships) on each simulated matrix with SBM algorithm.
#'
#' @return
#' A table with the memberships of each sample (in rows) as a function of each permuted matrix (in columns).
#' The number of columns corresponds to the number of permutations and the number of rows corresponds to the number of samples.
#'
#'
#' @export
Simulate_Clustering_SBM <- function(FB, NPERM, name) {

  plan(multisession, workers = parallel::detectCores() - 1)

  # Convert FB once (avoid making many copies in memory)
  FB_df <- as.data.frame(FB)

  # ---- STEP 1: parallel simulation ----
  Sim_one <- function(i, FB)
  {
    # One permutation only
    simulate(nullmodel(FB, "curveball"), nsim = 1)[,,1]
  }

  message("Running ", NPERM, " curveball permutations in parallel...")
  perm_list <- future_lapply(
    X = seq_len(NPERM),
    FUN = Sim_one,
    FB = FB_df
  )

  # perm_list is a list of NPERM matrices
  # Convert back to matrix 3D array if needed
  # OR keep list (better for parallel SBM)

  # ---- STEP 2: parallel SBM estimation ----
  SBM_one <- function(mat)
  {
    fit <- estimateBipartiteSBM(
      as.matrix(mat),
      model = 'bernoulli',
      estimOptions = list(verbosity = 0, plot = FALSE)
    )
    return(fit$memberships$row)
  }

  message("Running ", NPERM, " SBM fits in parallel...")
  cluster_list <- future_lapply(perm_list, SBM_one)

  # ---- STEP 3: Combine into dataframe ----
  simulated_clustering <- do.call(cbind, cluster_list)
  simulated_clustering <- as.data.frame(simulated_clustering)

  # ---- STEP 4: Export ----
  out_file <- file.path(
    here::here(),
    "data", "derived_data", "5_SBMclustering",
    name
  )

  write.table(
    simulated_clustering,
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  message("Saved simulated clustering to: ", out_file)

  return(simulated_clustering)
}
