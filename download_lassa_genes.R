library('tidyverse')
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library('genbankr')
library('ape')
library('rentrez')

set_entrez_key('2940c8fce95de8bb652be512772faf78a009')
##This is an API key that will allow you to download from genbank. You can set this up yourself or use mine.

accession_lassa <- read.csv('data/accession_lassa_all.seq', header = FALSE)
##All of the accession numbers that were captured in the search and downloaded to the accessionn_lassa_all.seq file have now been loaded in as accession_lassa
lassa_test_1 <- c('AF181853.1','KM821861.1','KM822082.1','MH887803.1',
                  'MK107932.1','MH053484.1','KM821847.1','AY628205.1')
##This is me just testing the process

accession_lassa <- GBAccession(accession_lassa$V1)
##We have now told R that they are accession numbers

Lassa_list <- list() #An empty list to be populated
for(entry in 1:length(accession_lassa)){
  entry.name <- accession_lassa[entry]
  Lassa_list[entry.name] <- readGenBank(GBAccession(accession_lassa[entry]), partial = T)
}

write_rds(Lassa_list, 'data/lassa_list.rds')
##So we don't need to download each time