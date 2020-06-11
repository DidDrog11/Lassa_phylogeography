library('tidyverse')
library('ggplot2')


lassa_list <- readRDS('data/lassa_list.rds')

##What happens here is that the code goes through each accession number in accession_lassa and tries to download it from genbank and puts it into the lassa_list
lassa_strains <- tibble(accession_lassa)
lassa_strains$host <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$host))

##Append where they were obtained from
lassa_strains$country <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$country))
lassa_strains$country <- gsub(':.*','\\',lassa_strains$country)

##Append when they were obtained
lassa_strains$year <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$collection_date))
RIGHT = function(x,n){
  substring(x,nchar(x)-n+1)
}
lassa_strains$year <- as.numeric(lassa_strains$year)

lassa_strain_x <- lassa_strains %>%
  filter(strain == x) %>%
  select(accession_numbers)

lassa_list[["AF182228.1"]]@sequence

lassa_strains$sequence <- for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]] %in% "lassa_strain_x")) {
    sequence <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]])
  myposition <- which(Lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "NP")
  myaa <- as.character(Lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
  Lassa_sequences$myaa[Lassa_sequences$accession_lassa_gb==entry] <- myaa
  Lassa_sequences$lmysea[Lassa_sequences$accession_lassa_gb==entry] <- nchar(myaa)
  }
}



lassa_mastomys <- lassa_strains %>%
  filter(host == 'Mastomys natalensis' | host == 'Mastomys')

lassa_hosts <- lassa_strains %>%
  filter(host != "NULL" & host != "Homo sapiens") %>%
  recode(lassa_hosts$host, 'Mastomys sp' = 'Mastomys sp.')

lassa_hosts$host <- lassa_hosts$host %>%
  recode(., 'Mastomys sp' = 'Mastomys sp.') %>%
  recode(., 'Mastomys' = 'Mastomys sp.')

#Plot the different species with detected lassa
species <- ggplot(data = lassa_hosts, aes(x = host))+
  geom_bar()+
  labs(title = "Non-human species providing Lassa sequences to GenBank",
       x = "Host Species", y = "Number of NCBI sequences")+
  theme_bw()

collected <- ggplot(data = lassa_strains, aes(x = year, fill = country))+
  geom_bar()+
  labs(title = "Year of Lassa sequences to GenBank",
       x = "Year", y = "Number of NCBI sequences")+
  theme_bw()
