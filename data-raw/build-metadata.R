# Run to update "survey metadata" from Excel file
# Saves under "data/survey_metadata.rda"

library(tidyverse)
library(readxl)
source('R/global_funcs.R')

survey_metadata <- read_xlsx(abs_pesp_path('Reference Documents/GroupMetadata.xlsx'), skip = 2)

usethis::use_data(survey_metadata, overwrite = TRUE)
