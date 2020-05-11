# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("genbankr")
library('genbankr')
# BiocManager::install('ape')
library('ape')
# install.packages('rentrez')
library('rentrez')
library('tidyverse')
library('ggplot2')

##Anything with install in it only needs to be run once. After it is installed you load it in using the library command. This needs to be done everytime you open R
##If you put a hash (#) at the beginning of a line it will prevent the line being run as code

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

##What happens here is that the code goes through each accession number in accession_lassa and tries to download it from genbank and puts it into the lassa_list
lassa_strains <- tibble(accession_lassa)
lassa_strains$host <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$host))

##Append where they were obtained from
lassa_strains$country <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$country))
lassa_strains$country <- gsub(':.*','\\',lassa_strains$country)

##Append when they were obtained
lassa_strains$year <- as.character(lapply(Lassa_list, function(x) x@sources@elementMetadata@listData$collection_date))
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}
lassa_strains$year <- RIGHT(lassa_strains$year, 4)


lassa_mastomys <- lassa_strains %>%
  filter(host == 'Mastomys natalensis' | host == 'Mastomys')
