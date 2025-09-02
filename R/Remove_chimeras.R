#' @title Identify and remove chimeras
#' @param fordada abundance table
#' @param name name will be used to import files and write names before extension of the different output tables
#' @description
#' - Tune data structure
#' - Remove chimeras with function isBimeraDenovo (package dada2)
#' - Export sequences that are chimeras ("_chimeras_list.tsv")
#' - Export summary file with number_of_total_clusters, number_of_bimeras, number_of_nobimeras ("_chimeras_summary.txt")
#'
#' @return
#' An abundance table without chimeras
#'
#' @examples
#'Remove_chimeras(fordada=fordada, name=name)
#'
#' @export
Remove_chimeras <- function(fordada, name){

  # Remove chimeras
  ## Tune data structure
  samples <- colnames(fordada)[seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))]
  print("samples:")
  print(samples)

  # Structuring the data table: seed_sequence in columns and samples in lines
  fordada.m <- data.matrix(fordada[ , seq(1+which(colnames(fordada)=="observation_sum"), ncol(fordada))])
  rownames(fordada.m) <- fordada$seed_sequence
  fordada.t <- t(fordada.m)

  ## Remove chimeras with function isBimeraDenovo
  bimeras.v <- dada2::isBimeraDenovo(dada2::getUniques(fordada.t), minFoldParentOverAbundance=4, verbose=TRUE)

  print(paste("Proportion of bimeras:", round(sum(as.logical(bimeras.v))/length(bimeras.v), 2)))

  bimeras.df <- data.frame(x=as.logical(bimeras.v), seq=names(bimeras.v))
  head(bimeras.df)

  ## Export sequences that are not chimeras
  # write.table(fordada[bimeras.df$x==FALSE,], file = file.path(here::here(), "data", "derived_data", paste(name,"_cleaned_abundance.tsv",sep="")), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)#writing of filtered abundance table (without identified chimeras)

  ## Export sequences that are chimeras
  write.table(fordada[bimeras.df$x==TRUE,], file = file.path(here::here(), "data", "derived_data", paste(name,"_chimeras_list.tsv", sep="")), sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)#downloading of identified chimeras

  ## Export summary file
  number_of_total_clusters = nrow(fordada)
  number_of_bimeras = nrow(fordada[bimeras.df$x==TRUE,])
  number_of_nobimeras = nrow(fordada[bimeras.df$x==FALSE,])
  summary.df = data.frame(number_of_total_clusters, number_of_bimeras, number_of_nobimeras)
  write.table(summary.df,file=file.path(here::here(), "data", "derived_data", paste(name,"_chimeras_summary.txt",sep="")),sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

  return(fordada[bimeras.df$x==FALSE,])
}
