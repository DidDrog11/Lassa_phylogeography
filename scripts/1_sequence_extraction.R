# The date we ran the search
search_date <- "2021-03-15"

# Set an entrez API key
set_entrez_key(entrez_key)

# We pull all the accession numbers found in a search for Lassa mammarenavirus on the nucleotide dataset from NCBI
ifelse(file.exists(here("data", "accession_2021-03-15.rds")),
        accession_lassa <- read_rds(here("data", "accession_2021-03-15.rds")),
        accession_lassa <- searchGB(query = '"Lassa mammarenavirus"[Organism]', db = "nucleotide", sequences = F)
        )

# We save the accession numbers we have identified on the day we ran the search
write_rds(accession_lassa, file = here("data", paste("accession_", search_date, ".rds", sep = "")))

# We define these as accession numbers
accession_lassa <- GBAccession(accession_lassa)

# If we already have the data from GenBank we read it in here, otherwise we download it
ifelse(file.exists(here("data", "lassa_list_2021-03-15.rds")),
       lassa_list <- read_rds(here("data", "lassa_list_2021-03-15.rds")),
       lassa_list <- download_lassa(accession_lassa) # This uses the readGenBank function to download the entries for each accession number in accession_lassa and save the output
       )

# This produces a table of all the accession numbers that we will fill with information we want
lassa_strains <- tibble(accession_lassa)
# This labels the accession number based on whether lassa is mentioned in the entries definition
lassa_strains$defined_lassa <- as.character(lapply(lassa_list, function(x) x@definition)) %>%
  str_detect(., "(Lassa)|(lassa)")
# We extract the host species
lassa_strains$host <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$host))
# Append where they were obtained from at country and region level
lassa_strains <- tibble(as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$country))) %>%
  separate(., col = 1, into = c("country", "region"), sep = ": ") %>%
  cbind(lassa_strains, .)
# Append when they were obtained
lassa_strains$year <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$collection_date)) %>%
  lubridate::parse_date_time(., c("%d-%b-%Y", "%b-%Y", "%Y")) %>%
  lubridate::year(.)
# We use the lubridate package to convert the dates into the same format and then extract the year for each entry

# The names of the genes saved on GenBank
lassa_strains$gene <- as.character(lapply(lassa_list, function(x) x@cds@elementMetadata@listData$gene))
unique(lassa_strains$gene) # These sequences have been uploaded as GPC, NP, Z, L, pol etc. The columns containing a c means the accession number contains sequences for 2 genes
lassa_strains$notes <- as.character(lapply(lassa_list, function(x) x@definition))

# We extract the nucleotide sequences from the genbank entries for the NP, GPC, Z and L 

## This extracts all the nucleoprotein sequences we create two variables the sequence and the length
all_sequences <- lassa_strains
all_sequences$np_seq <- NA
all_sequences$np_length <- NA

for(entry in names(lassa_list)) {
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("NP", "np", "N"))) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("NP", "np", "N"))
    np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
    all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
  } else { # These first three statements look in the gene section of the record if it is NP, np or N it is extracted
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% "nucleoprotein")
      np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
      all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
    } else { # This statement looks in the note, if it includes nucleoprotein it is extracted
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$product %in% "nucleoprotein")
        np_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
        all_sequences$np_seq[all_sequences$accession_lassa==entry] <- np_seq
        all_sequences$np_length[all_sequences$accession_lassa==entry] <- nchar(np_seq)
      } else { # This looks within the product section
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
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("GP", "GPC", "gpc"))) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("GP", "GPC", "gpc"))
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
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("L", "l"))) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("L", "l"))
    l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
    all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% c("L", "l"))) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% c("L", "l"))
      l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
      all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
    } else {
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "l protein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "l protein")
        l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
        all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
        all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "L protein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "L protein")
          l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
          all_sequences$l_seq[all_sequences$accession_lassa==entry] <- l_seq
          all_sequences$l_length[all_sequences$accession_lassa==entry] <- nchar(l_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "polymerase")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% c("L polymerase", "polymerase"))
            l_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
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
  if(any(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("Z", "z"))) {
    myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$gene %in% c("Z", "z"))
    z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
    all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
    all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
  } else {
    if(any(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% c("Z", "z"))) {
      myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData$note %in% c("Z", "z"))
      z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData$translation[myposition])
      all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
      all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
    } else {
      if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "z protein")) {
        myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "z protein")
        z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
        all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
        all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
      } else {
        if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "Z protein")) {
          myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "Z protein")
          z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
          all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
          all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
        } else {
          if(any(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "ring-finger protein")) {
            myposition <- which(lassa_list[[entry]]@cds@elementMetadata@listData[["product"]] %in% "ring-finger protein")
            z_seq <- as.character(lassa_list[[entry]]@cds@elementMetadata@listData[["translation"]][myposition])
            all_sequences$z_seq[all_sequences$accession_lassa==entry] <- z_seq
            all_sequences$z_length[all_sequences$accession_lassa==entry] <- nchar(z_seq)
          }
        }
      }
    }
  }
}

# We then extract the nucleotide sequence from entry once we know what it represents
all_sequences$name_nuc_seq <- NA
all_sequences$nuc_seq <- NA
all_sequences$nuc_length <- NA

for(entry in names(lassa_list)) {
    name_nuc <- as.character(names(lassa_list[[entry]]@sequence))
    nuc_seq <- as.character(unname(lassa_list[[entry]]@sequence))
    all_sequences$name_nuc_seq[all_sequences$accession_lassa==entry] <- name_nuc
    all_sequences$nuc_seq[all_sequences$accession_lassa==entry] <- nuc_seq
    all_sequences$nuc_length[all_sequences$accession_lassa==entry] <- nchar(nuc_seq)
    }

all_sequences %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  select(np_length, gpc_length, l_length, z_length) %>%
  rowSums(.) -> all_sequences$all_seqs

# These accession numbers are genbank sequences where we have been unable to extract data for the GPC, NP, L or Z
missing_sequences <- all_sequences %>%
  filter(all_sequences$all_seqs == 0 & all_sequences$defined_lassa == T) %>%
  select(accession_lassa, gene, notes)
missing_seq <- as.character(missing_sequences$accession_lassa)
missing_sequences_list <- lassa_list[c(missing_seq)]

write_rds(all_sequences, here("data", paste("all_sequences_", search_date, ".rds", sep = "")))