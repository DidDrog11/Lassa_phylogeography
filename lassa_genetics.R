library('tidyverse')
library('genbankr')
library('ggplot2')


lassa_list <- readRDS('data/lassa_list.rds')
accession_lassa <- read.csv('data/accession_lassa_all.seq', header = FALSE)

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

##What data do we have on the strain from GenBank?
lassa_strains$gene <- as.character(lapply(lassa_list, function(x) x@cds@elementMetadata@listData$gene))
unique(lassa_strains$gene) ##We have 20 different "strains" according to how people entered them. The ones with c infront means we have two sequences for that accession number
lassa_strains$notes <- as.character(lapply(lassa_list, function(x) x@definition))

## This extracts all the nucleoprotein sequences
all_sequences <- lassa_strains %>%
  rename("accession_lassa" = 1)
all_sequences$np_seq <- NA
all_sequences$np_length <- NA

for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "NP")) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "NP")
    np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
    all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "np")) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "np")
      np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
      all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
    } else {
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")
        np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
        all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
        all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")
          np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
          all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
          all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "nucleocapsid protein")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "nucleocapsid protein")
            np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
            all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
            all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
          } else {
            if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "hypothetical protein")) {
              myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "hypothetical protein")
              np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
              all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
              all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
            } else { 
              np_seq <- as.character(lassa_list[["X52400.1"]]@cds@elementMetadata@listData$translation[2])
              all_sequences$np_seq[all_sequences$accession_lassa=="X52400.1"] <- np_seq
              all_sequences$np_length[all_sequences$accession_lassa=="X52400.1"] <- nchar(np_seq)
            }
          }
        }
      }
    }
  }
}

all_sequences$gpc_seq <- NA
all_sequences$gpc_length <- NA
## This extracts all the glycoprotein sequences
for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "GP")) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "GP")
    gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
    all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "GPC")) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "GPC")
      gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
      all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
    } else {
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "gpc")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "gpc")
        gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
        all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
        all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "glycoprotein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "glycoprotein")
          gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
          all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
          all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "glycoprotein precursor")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "glycoprotein precursor")
            gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
            all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
            all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
          } else {
            if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "glycoprotein")) {
              myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "glycoprotein")
              gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
              all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
              all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
            } else {
              if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "glycoprotein precursor")) {
                myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "glycoprotein precursor")
                gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
                all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
                all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
              } else {
                if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "envelope glycoprotein")) {
                  myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "envelope glycoprotein")
                  gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
                  all_sequences$gpc_seq[all_sequences$accession_lassa==entry] <- gpc_seq
                  all_sequences$gpc_length[all_sequences$accession_lassa==entry] <- nchar(gpc_seq)
                } else { 
                  gpc_seq <- as.character(lassa_list[["X52400.1"]]@cds@elementMetadata@listData$translation[1])
                  all_sequences$gpc_seq[all_sequences$accession_lassa=="X52400.1"] <- gpc_seq
                  all_sequences$gpc_length[all_sequences$accession_lassa=="X52400.1"] <- nchar(gpc_seq)
                }
              }
            }
          }
        }
      }
    }
  }
}

all_sequences$l_seq <- NA
all_sequences$l_length <- NA
## This extracts all the L
for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "L")) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "L")
    l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
    all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "L")) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "L")
      l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
      all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
    } else{
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "l protein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "l protein")
        gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
        all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
        all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "L protein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "L protein")
          gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
          all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
          all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "polymerase")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "polymerase")
            gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
            all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
            all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
          }
        }
      }
    }
  }
}

all_sequences$z_seq <- NA
all_sequences$z_length <- NA
## This extracts all the Z
for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "Z")) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "Z")
    z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
    all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "z")) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% "z")
      z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
      all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
    } else {
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "z protein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "z protein")
        gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
        all_sequences$z_seq[all_sequences$accession_lassa==entry] <- l_seq
        all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "Z protein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "Z protein")
          gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
          all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
          all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "ring-finger protein")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "ring-finger protein")
            gpc_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
            all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
            all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
          }
        }
      }
    }
  }
}


all_sequences %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(np_length, gpc_length, l_length, z_length) %>%
  rowSums(.) -> all_sequences$all_seqs

missing_sequences <- all_sequences %>%
  filter(all_sequences$all_seqs == 0) %>%
  select(accession_lassa, gene, notes)
missing_seq <- as.character(missing_sequences$accession_lassa)
missing_sequences_list <- lassa_list[c(missing_seq)]

all_sequences[is.na(all_sequences)] <- 0

all_seq <- all_sequences
all_sequences <-  replace(is.na(select_if(all_sequences, is.numeric)), 0)
  replace(is.na(.), 0) %>%
  mutate(all_aa_lengths = rowSums(.))
  mutate(all_aa_lengths = sum(
    replace(
      is.na(select(np_length, gpc_length, l_length, z_length)
            , 0))))
  replace(is.na(.), 0) %>%
  mutate(all_aa_lengths = sum(p_length + gpc_length + l_length + z_length))
  mutate(all_aa_lengths = rowSums(np_length + gpc_length + l_length + z_length, na.rm = T))

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
