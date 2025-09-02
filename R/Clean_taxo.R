#' @title Clean the taxonomy
#' @param table multi-affiliation or abundance table with "bidouille" column (a copy of the blast_taxonomy column)
#' @description
#'  Clean the taxonomy (7 ranks: Kingdom, Phylum, Class, Order, Family, Genus, Species) multi-affiliation or abundance table and follow these rules:
#'  - Delete a lowercase followed by two underscores
#'  - Replace an underscore with a space
#'  - Delete for the special characters []()_"0123456789.
#'  - Replace all occurrences of the "one or more spaces followed by a semicolon" pattern with a single semicolon
#'  - Replacement of 'Chloroplast' and 'Mitochondria' by 'Plantae_Chl' and 'Plantae_Mit' respectively in all 7 taxonomic fields
#'  - Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields
#'  - Overwriting of the genus field with 'unknown genus'  if presence of 'incertae sedis' in any of upper taxonomic fields
#'  - Remove words after the first word found in the 6 first taxonomic fields (except if 'unknown' is found)
#'  - Remove words after the second word found in the 7th taxonomic field
#'  - Replacement in the species field of 'unidentified' or 'metagenome' or 'endosymbiont' or 'organism' anywhere in the field by 'unknown species'
#'  - Replacement in the species field of 'sp', 'eae bacterium', 'ales bacterium', 'degrading'  at the end of the field by 'unknown species'
#'  - Replacement in the species field of 'bacterium' or 'uncultured' at the beginning of the field by 'unknown species'
#'  - Replacement of the species field by 'unknown species' if empty
#'  - Overwriting of the species field with 'unknown species' if it begins by a lower case or contains " bacterium"
#'  - Overwriting of the genus field  (if unknown) with the first word of the species field (if not unknown) if empty
#'
#' @return
#' An table with a clean taxonomy : 7 columns (tx1, tx2, tx3, tx4, tx5, tx6, tx7) for the 7 taxonomic ranks, replacing the "bidouille" column
#'
#' @examples
#' my_ab_result <- Clean_taxo(table = ab)
#'
#' @export
Clean_taxo <- function(table)
{

  ## Delete a lowercase followed by two underscores
  table$bidouille <- gsub("[a-z]__", "", table$bidouille, perl=TRUE)
  ## Replace an underscore with a space
  table$bidouille <- gsub("_", " ", table$bidouille, perl=TRUE)
  ## Delete for the special characters []()_"0123456789.
  table$bidouille <- gsub("[\\[\\]()_\"0-9\\.]", "", table$bidouille, perl=TRUE)
  ## Replace all occurrences of the "one or more spaces followed by a semicolon" pattern with a single semicolon
  table$bidouille <- gsub(" +;", ";", table$bidouille, perl=TRUE)

  ## Deconstruction of the taxonomic column ('blast_taxonomy') in 7 taxonomic fields
  ## Warnings message when no data on taxonomy
  table[table$bidouille=="no data",]
  table <- table %>% tidyr::separate(col=bidouille, into=c("tx1","tx2", "tx3", "tx4", "tx5", "tx6", "tx7"), sep=";")

  ## Replacement of 'Chloroplast' and 'Mitochondria' by 'Plantae_Chl' and 'Plantae_Mit' respectively in all 7 taxonomic fields
  table[grepl("Chloroplast", table$tx4, ignore.case=TRUE), 1:7] <- c("Plantae_Chl")
  table[grepl("Mitochondria", table$tx5, ignore.case=TRUE), 1:7] <- c("Plantae_Mit")

  ## Overwriting of the species field with 'unknown species' if presence of 'incertae sedis' in any of upper taxonomic fields
  ## Overwriting of the genus field with 'unknown genus' if presence of 'incertae sedis' in any of upper taxonomic fields
  ## Overwriting in the first 5 ranks 'unknown' if presence of 'incertae sedis'
  table[!(grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx2, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx3,ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx4, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx5, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx6, ignore.case=TRUE))&grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx7, ignore.case=TRUE), 7]<-c("unknown species")
  table[grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx2, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx3, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx4, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx5, ignore.case=TRUE)|grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", table$tx6, ignore.case=TRUE), 6:7]<-list("unknown genus","unknown species")
  table[, 1:5] <- lapply(table[, 1:5], function(x) ifelse(grepl("incertae sedis|Incertae sedis|Incertaesedis|Incertae Sedis", x, ignore.case = TRUE), "unknown", x))

  ## Remove words after the first word found in the 6 first taxonomic fields (except if 'unknown' is found)
  ## Removal of 'Candidatus' (i.e.non cultivable) in each field
  table <- as.data.frame(table)
  for(j in 1:6) {
    table[, j] <- gsub("^Candidatus ", "", table[,j ])
    my.stuff <- as.data.frame(table[, j])
    tf <- !grepl("unknown", my.stuff[, 1], ignore.case=TRUE)
    table[which(tf),j] <- gsub(" +.*$", "", table[which(tf), j], perl=TRUE)
  }
  table <-  data.table::as.data.table(table)

  ## Cleaning for more than 2 words in the 7th taxonomic field
  ## Removal of 'Candidatus' (i.e.non cultivable) in each field
  table$tx7 <- gsub("^Candidatus ","", table$tx7, perl=TRUE)
  ## Remove words after the second word found in the 7th taxonomic field
  table$tx7 <- gsub(" +$", "", gsub("^(([^ ]+ ?){1,2}).*$","\\1",table$tx7,perl=TRUE),perl=TRUE)

  ## Replacement in the species field of 'unidentified' or 'metagenome' or 'endosymbiont' or 'organism' anywhere in the field by 'unknown species'
  table$tx7 <- gsub("^.*(unidentified|metagenome|endosymbiont|organism).*$","unknown species",table$tx7,perl=TRUE)
  ## Replacement in the species field of 'sp', 'eae bacterium', 'ales bacterium', 'degrading'  at the end of the field by 'unknown species'
  table$tx7 <- gsub("^.*(sp|eae bacterium|ales bacterium|degrading)$","unknown species",table$tx7,perl=TRUE)
  ## Replacement in the species field of 'bacterium' or 'uncultured' at the beginning of the field by 'unknown species'
  table$tx7 <- gsub("^(uncultured|bacterium).*$","unknown species", table$tx7, perl=TRUE)
  ## Replacement of the species field by 'unknown species' if empty
  table$tx7[grepl("^ *$", table$tx7)] <-"unknown species"

  ## Overwriting of the species field with 'unknown species' if it begins by a lower case or contains " bacterium"
  todo <- grepl("^[a-z]| bacterium", table$tx7) & !grepl("unknown species", table$tx7)
  if(sum(todo)!=0){
    print("Replacing the following species names by \'unknown species\' ")
    print( table$tx7[which(todo)])
  }
  table$tx7[which(todo)] <- "unknown species"

  ## Overwriting of the genus field  (if unknown) with the first word of the species field (if not unknown) if empty
  todo <- table$tx6=="unknown genus" & grepl(" ", table$tx7) & !grepl("unknown species", table$tx7)
  if(sum(todo)!=0){
    print("Replacing the following genus names ...")
    print( table$tx6[which(todo)])
    print("... with the first word of the following species name")
    print( table$tx7[which(todo)])
  }
  table$tx6[which(todo)] <- gsub(" .*$", "", table$tx7[which(todo)])

  return(table)

}

