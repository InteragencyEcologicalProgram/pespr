# Run to update "phyto taxonomy" from csv data
# Saves under "data/phyto_taxonomy.rda"

library(tidyverse)
source('R/global_funcs.R')

csv_path <- abs_pesp_path('Reference Documents/PhytoTaxonomy.csv')

phyto_taxonomy <- read_quiet_csv(csv_path)

usethis::use_data(phyto_taxonomy, overwrite = TRUE)

attr(phyto_taxonomy, 'version')     <- '2026.1.1'
attr(phyto_taxonomy, 'updated')     <- as.character(Sys.Date())