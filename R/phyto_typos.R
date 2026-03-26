#' Phytoplankton taxon name typo lookup table
#'
#' A reference table of known misspellings and non-standard taxon name variants
#' used to correct raw data prior to synonym resolution. Created by DWR staff.
#' Used internally by `correct_taxon_typos()`.
#'
#' @format A dataframe with columns:
#' \describe{
#'   \item{Taxon}{The misspelled or non-standard taxon name as it appears in raw data}
#'   \item{TaxonCorrected}{The corrected taxon name to substitute}
#' }
"phyto_typos"
