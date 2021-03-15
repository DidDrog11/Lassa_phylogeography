download_lassa <- function(accession_numbers){
  lassa_list <- list() # An empty list to be populated
  # For each accession number we download the relevant entry from genbank
  for(entry in 1:length(accession_lassa)){
    entry.name <- accession_lassa[entry]
    lassa_list[entry.name] <- readGenBank(GBAccession(accession_lassa[entry]), partial = T)
  } # There are over 2,000 sequences so this takes a while
  write_rds(lassa_list, here("data", paste("lassa_list_", search_date, ".rds", sep = ""))) # Once downloaded the output is saved so we don't need to download again
}

