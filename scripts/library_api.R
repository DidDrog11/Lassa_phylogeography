# We load packages we need
if (!require("pacman")) install.packages("pacman")
pkgs = 
  c("here",
    "cowplot",
    "conflicted",
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
    "readxl",
    "mgcv",
    "mgcViz",
    "raster",
    "colorspace",
    "scatterpie"
  )
pacman::p_load(pkgs, character.only = T)

devtools::install_github("yutannihilation/ggsflabel")
library("ggsflabel")

# We use an API to download the data from LASV related accession codes
entrez_key <- rstudioapi::askForSecret("entrez API key")
google_key <- rstudioapi::askForSecret("google API key")

register_google(key = google_key)

conflict_prefer("select", "dplyr", "raster")
conflict_prefer("filter", "dplyr")
conflict_prefer("getData", "raster")
