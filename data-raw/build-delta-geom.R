# data-raw/build-delta-geom.R
# purpose: build delta_geom dataset for pespr
# source: deltamapr::R_EDSM_Subregions_Mahardja
# retrieved: 2026-01-23
# deltamapr version: 1.0.1
# notes: regional boundary designations for the Delta (Mahardja subregions)

library(deltamapr)

delta_geom <- deltamapr::R_EDSM_Subregions_Mahardja

attr(delta_geom, 'source') <- 'deltamapr::R_EDSM_Subregions_Mahardja'
attr(delta_geom, 'source_version') <- as.character(utils::packageVersion('deltamapr'))
attr(delta_geom, 'retrieved') <- as.character(Sys.Date())

usethis::use_data(delta_geom, overwrite = TRUE)

usethis::use_data_doc('delta_geom')
