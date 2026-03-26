#' Survey program metadata
#'
#' Metadata for PESP survey programs.
#' Used internally by `read_meta_file()` and `add_meta_col()`.
#'
#' @format A dataframe with columns including:
#' \describe{
#'   \item{Survey}{Program identifier (e.g., `'EMP'`)}
#'   \item{Starting Date}{Start date for the metadata record}
#'   \item{Ending Date}{End date for the metadata record (`NA` treated as current date at runtime)}
#' }
"survey_metadata"
