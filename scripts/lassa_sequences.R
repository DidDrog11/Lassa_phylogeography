library('tidyverse')
library('ggplot2')
library('ggmap')
library('sf')
library("rnaturalearth")
library("rnaturalearthdata")


lassa_list <- readRDS('data/lassa_list.rds')
accession_lassa <- read.csv('data/accession_lassa_all.seq', header = FALSE)

##What happens here is that the code goes through each accession number in accession_lassa and tries to download it from genbank and puts it into the lassa_list
lassa_strains <- tibble(accession_lassa)
lassa_strains$host <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$host))

##Append where they were obtained from
lassa_strains$country <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$country))
lassa_strains$country <- gsub(':.*','\\',lassa_strains$country)
lassa_strains$region <- as.character(lapply(lassa_list, function(x) x@sources@elementMetadata@listData$country)) %>%
  gsub('.*: ','',.)
lassa_strains <- lassa_strains %>%
  mutate(region = ifelse(region %in% country, NA, region),
         region = ifelse(region == "NULL", NA, region),
         country = ifelse(country == "NULL", NA, country),
         location = ifelse(is.na(region), NA, paste(region, country, sep = ", ")))

geocoded_region <- na.omit(unique(lassa_strains$location)) %>%
  tibble() %>%
  rename(location = 1) %>%
  mutate_geocode(location)

# This misses some places for example Bantou, Guinea we can extract these from searching on google maps
# Jirandogo and Tanganya are estimated
missing_coords <- tibble(location = c("Bantou, Guinea", "Brissa, Guinea", "Denguedou, Guinea", "Gbetaya, Guinea", "Jirandogo, Ghana", "Khoria, Guinea", "Tanganya, Guinea", "Worogui, Benin",
                                      "Safrani, Guinea", "Silimi, Guinea", "Odo-akaba, Benin"),
                         lon = c(10.05, 10.22, 8.49, 9.84, 8.65, 9.94, 9.88, 8.88, 0.06, 9.98, 8.77),
                         lat = c(-10.59, -10.69, -10.45, -11.04, 0.00, -10.89, -10.8, -2.67, -10.74, -10.82, -2.6))

locations <- geocoded_region %>%
  left_join(., missing_coords, by = "location") %>%
  mutate(lon = coalesce(lon.x, lon.y),
         lat = coalesce(lat.x, lat.y)) %>%
  select(location, lon, lat)

geocoded_nation <- na.omit(unique(lassa_strains$country)) %>%
  tibble() %>%
  rename(country = 1) %>%
  mutate_geocode(country) %>%
  mutate(lon = ifelse(country == "Togo", 8.6195, lon),
         lat = ifelse(country == "Togo", 0.8248, lat))

lassa_strains <- lassa_strains %>%
  left_join(., locations, by = "location") %>%
  left_join(., geocoded_nation, by = "country") %>%
  mutate(lon = ifelse(is.na(lon.x), lon.y, lon.x),
         lat = ifelse(is.na(lat.x), lat.y, lat.x)) %>%
  select(V1, host, country, region, location, lon, lat)

boundaries <- tibble(lon_max = max(na.omit(lassa_strains$lon))+1,
                     lon_min = min(na.omit(lassa_strains$lon))-1,
                     lat_max = max(na.omit(lassa_strains$lat))+1,
                     lat_min = min(na.omit(lassa_strains$lat))+1)

world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(boundaries$lon_min, boundaries$lon_max), ylim = c(boundaries$lat_min, boundaries$lat_max), expand = FALSE) +
  geom_count(data = lassa_strains, aes(x = lon, y = lat)) +
  theme_bw()

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