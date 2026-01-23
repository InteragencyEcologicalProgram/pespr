# Run to update "taxa typos" from csv data
# Saves under "data/phyto_typos.rda"

library(tidyverse)
source('R/global_funcs.R')

csv_path <- abs_pesp_path('Reference Documents/TaxaTypos.csv')

phyto_typos <- read_quiet_csv(csv_path)

usethis::use_data(phyto_typos, overwrite = TRUE)