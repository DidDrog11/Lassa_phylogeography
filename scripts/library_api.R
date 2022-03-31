# We load packages we need
if (!require("pacman")) install.packages("pacman")
pkgs = 
  c("here",
    "cowplot",
    "genbankr", 
    "ape",
    "rentrez",
    "tidyverse",
    "ggplot2",
    "ggspatial",
    "insect",
    "lubridate",
    "tmap",
    "ggmap",
    "sf",
    "rnaturalearth",
    "rnaturalearthdata",
    "readxl"
  )
pacman::p_load(pkgs, character.only = T)

# We use an API to download the data from LASV related accession codes
entrez_key <- rstudioapi::askForSecret("entrez API key")
