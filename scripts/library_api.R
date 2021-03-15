# Biocmanager handles the installation and updating of some of the packages we use
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

BiocManager::install(pkgs = c("genbrankr", "ape", "rentrez"))

# We install other packages of interest
if (!require("pacman")) install.packages("pacman")
pkgs = 
  c("genbankr", 
    "ape",
    "rentrez",
    "tidyverse",
    "ggplot2",
    "insect",
    "lubridate"
  )
pacman::p_load(pkgs, character.only = T)

# We use an API to download the data from LASV related accession codes
rstudioapi::askForSecret("entrez API key")