#' Phytoplankton taxonomy reference table
#'
#' Hierarchical taxonomic classifications for phytoplankton taxa observed in
#' PESP programs based on AlgaeBase. Algal group assignments from Tiffany
#' Brown (DWR). Includes synonym mapping via the `CurrentTaxon` column.
#'
#' @format A dataframe with columns:
#' \describe{
#'   \item{Kingdom}{Kingdom-level classification}
#'   \item{Phylum}{Phylum-level classification}
#'   \item{Class}{Class-level classification}
#'   \item{AlgalGroup}{Functional algal group (DWR classification)}
#'   \item{Genus}{Genus name}
#'   \item{Species}{Species}
#'   \item{Taxon}{Accepted taxon name}
#'   \item{CurrentTaxon}{Current accepted name if `Taxon` is a synonym; `'None'` if already current}
#' }
#' @source \url{https://www.algaebase.org}
"phyto_taxonomy"
