#' @title Define multi-affiliation
#' @param multiaffiliation_table multi-affiliation table obtained from FROGS
#' @param thresh frequency threshold below which 'unknown genus' or 'unknown species' are not considered to infer multi-affilitation
#' @description
#' Work on cleaning up the taxonomy and define a taxonomy affiliation from clusters/OTUs with several affiliations
#'- Cleaning up the taxonomy using function Clean_taxo
#'
#'- Define a taxonomy affiliation from clusters/OTUs with several affiliations :
#'define a threshold (tresh) for which, in each cluster, affiliations with "unknown genus" or "unknown species" and their proportion is lower than the threshold, are not considered to define the affiliation of the cluster/OTU (column ignore_genus and ignore_species)
#'
#' For a given taxonomic level :
#'- consider only the names for the species names which are not ignored
#'- consider all the names, only if all species names are ignored
#'- if only one taxonomic name, keep it
#'- if two or more taxonomic names, keep unknown name (if it present) or Multi-affiliation
#'
#'
#' @return
#' An multi-affiliation table with a unique taxonomy affiliation per clusters/OTUs
#'
#' @examples
#' taxo_16S_Multiaffiliation <- Define_multiaffiliation(multiaffiliation_table = ma_MS_16S, thresh = 0.5)
#'
#' @export
Define_multiaffiliation <- function(multiaffiliation_table, thresh)
{

  # Import Multi-affiliation data
  mt <-  multiaffiliation_table
  names(mt)
  # Create bidouille column (a copy of the blast_taxonomy column)
  mh <- data.table::data.table(bidouille=mt$blast_taxonomy, blast_taxonomy=mt$blast_taxonomy, observation_name=mt$`#observation_name`)
  names(mh)

  # Remove file
  rm(mt)

  # Clean the taxonomy in mh (on bidouille column)
  mhh <- Clean_taxo(table =mh)

  # Ignoring of 'unknown species' or 'unknown genus' under a threshold
  ## Add a column with the number of the cluster (clus) and remove blast_taxonomy column
  mhh$clus <- as.integer(gsub("Cluster_", "", mhh$observation_name))
  mhh$blast_taxonomy <- NULL

  #### Add a column with the number of rows with the same affiliation until genus by cluster (n) and the total of affiliation by cluster (ntot) in multi-affiliation table
  ## For each cluster, count the number of rows with the exact same taxonomy until the genus level
  my_genus_result <- data.table::as.data.table(unique(mhh %>% dplyr::add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, observation_name)))
  ## For each cluster, count the total number of rows
  my_genus_result$ntot <-  data.table::as.data.table(unique(mhh %>% dplyr::add_count(observation_name)))$n
  ## Sort data.table
  data.table::setkey(my_genus_result, clus, n)
  my_genus_result$prop_genus <- my_genus_result$n/my_genus_result$ntot
  ## If genus is unknown and below a threshold, it is ingored
  my_genus_result$ignore_genus <- FALSE
  my_genus_result$ignore_genus <- my_genus_result$ignore_genus | (grepl("unknown", my_genus_result$tx6, perl=TRUE, ignore.case=TRUE) & my_genus_result$prop_genus < thresh)

  ## Add a column with the number of rows with the same affiliation until species by cluster (n) and the total of affiliation by cluster (ntot) in multi-affiliation table
  ## For each cluster, count the number of rows with the exact same taxonomy until the species level
  my_species_result <- data.table::as.data.table(unique(mhh %>% dplyr::add_count(clus, tx1, tx2, tx3, tx4, tx5, tx6, tx7, observation_name)))
  my_species_result$ntot <- data.table::as.data.table(unique(mhh %>% dplyr::add_count(observation_name)))$n
  ## Sort data.table
  data.table::setkey(my_species_result, clus, n)
  my_species_result$prop_species <- my_species_result$n/my_species_result$ntot
  ## If species is unknown and below a threshold, it is ingored
  my_species_result$ignore_species <- FALSE
  my_species_result$ignore_species <- my_species_result$ignore_species | (grepl("unknown", my_species_result$tx7, perl=TRUE,ignore.case=TRUE) & my_species_result$prop_species < thresh)

  # Creates a dataset with the taxonomy, cluster name and the columns "ignore_genus" and "ignore_species"
  my.first.result <- merge(my_genus_result, my_species_result, by = names(my_species_result)[1:8])[, c(1:8,13,18)]

  ## Check that no genus ignored when species is not ignored
  # sum(my.first.result[my.first.result$ignore_genus != my.first.result$ignore_species]$ignore_genus)

  ## Check that no species ignored when genus is not ignored
  # sum(my.first.result[my.first.result$ignore_species != my.first.result$ignore_genus]$ignore_species)

  # example : my.first.result <- my.first.result[my.first.result$observation_name == "Cluster_174",]

  ## Create a dataset with one row per cluster
  my.second.result <- data.table::data.table(observation_name=unique(my.first.result$observation_name))

  ## Loop over each taxonomic levels above the genus level
  for(j in 1:7) {
    #for(j in 1:5) {
    tx <- as.character(as.data.frame(my.first.result)[,j])
    print(tx)
    tmp1 <- apply(my.second.result, 1, function(r){
      if(j > 5){
        if(j==6){
          v <- tx[ !my.first.result$ignore_genus & my.first.result$observation_name == r[1]]
        }else{
          v <- tx[ !my.first.result$ignore_species & my.first.result$observation_name == r[1]]
        }
        #v <- ifelse(j==6, tx[ !my.first.result$ignore_genus & my.first.result$observation_name=="Cluster_247"], tx[ !my.first.result$ignore_species & my.first.result$observation_name==r[1]])
        #  unique(tx[(j==6 | !my.first.result$ignore_genus) & my.first.result$observation_name==r[1]])

        #v <- unique(tx[(j==7 | !my.first.result$ignore_species) & my.first.result$observation_name==r[1]])

      }else{

        ## For a given taxonomic level, consider all the names, only if all species names are ignored
        if(sum(my.first.result$ignore_species[my.first.result$observation_name==r[1]])==sum(my.first.result$observation_name==r[1])){
          v <- unique(tx[my.first.result$observation_name==r[1]])
        }else{
          ## otherwise, consider only the names for the species names which are not ignored
          v <- unique(tx[!my.first.result$ignore_species & my.first.result$observation_name==r[1]])
        }

      }
      print(v)
      ## If only one taxonomic name keeps it, otherwise "Multi-affiliation"
      if (length(unique(v))==1) {v[1]}
      else if ((length(unique(v))>=2)) {
        ## If two taxonomic names keeps unkown or Multi-affiliation
        unk <- grepl("unknown",v)
        unkV <- unique(v[unk])
        if (length(unkV) == 0) {
          "Multi-affiliation"
        } else {
          unkV[1]
        }
      }

      else {
        "Multi-affiliation"
      }

    })
    my.second.result <- cbind(my.second.result, tmp1)
  }

  # Reconstruction of the new full taxonomy from each taxonomic field for each cluster
  my_mhh_result <- data.table::data.table(observation_name=my.second.result$observation_name,  blast_taxonomy= apply(my.second.result, 1, function(r) {
    paste(r[2:8], sep=";", collapse=";")
  }))

  return(my_mhh_result)

}
