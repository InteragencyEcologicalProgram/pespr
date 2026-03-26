
# Read in Data ------------------------------------------------

# Reads a BSA Excel file, parsing Date/Time columns as date type
#' @importFrom readxl read_excel
#' @noRd
read_bsa_xlsx <- function(path){
  # read in datetime columns as date
  vec_names <- names(read_excel(path, n_max = 0))

  # define col types
  vec_types <- ifelse(grepl('Date|Time', vec_names), 'date', 'guess')

  # read in data
  df <- readxl::read_excel(path, col_types = vec_types)

  return(df)
}

# Reads a BSA CSV file, parsing Date/Time columns as date type
#' @importFrom readr read_csv cols col_date
#' @noRd
read_bsa_csv <- function(path){
  # read in column names to identify date/time columns
  vec_names <- names(read_csv(path, n_max = 0, show_col_types = FALSE))

  # define col types
  date_cols <- vec_names[grepl('Date|Time', vec_names)]
  col_spec <- do.call(cols, setNames(
    lapply(date_cols, function(x) col_date()),
    date_cols
  ))

  # read in data
  df <- read_csv(path, col_types = col_spec)

  return(df)
}

# Standardize Column Names ------------------------------------------------

# Helper Functions --------------------------------------------------------

# Drops Dimension and Measurement columns and removes blank rows
#' @importFrom dplyr select
#' @noRd
remove_rows_bsa <- function(df){
  # subset columns with word 'measurement'
  mes_cols <- colnames(df)[grepl('Measurement', colnames(df))]

  # remove measurement and dimension columns
  df <- df %>%
    select(-c(Dimension, mes_cols))

  # remove completely blank rows cause by measurement columns
  df <- janitor::remove_empty(df, which = 'rows')

  return(df)
}

# Renames BSA-specific col names to PESP col names and adds Lab = 'BSA'
#' @importFrom dplyr rename mutate
#' @noRd
rename_cols_bsa <- function(df){
  df <- df %>%
    rename(
      Date = SampleDate,
      Time = SampleTime,
      Station = StationCode,
      PhytoForm = `Colony/Filament/Individual Group Code`,
      GALD = `GALD 1`,
    ) %>%
    mutate(
      Lab = 'BSA'
    )

  return(df)
}

#' @title Calculate per-mL density values from BSA data
#'
#' @description
#' Computes standardized density columns (`Units_per_mL`, `Cells_per_mL`,
#' and/or `Biovolume_per_mL`) by multiplying raw counts by a factor column.
#' Biovolume is calculated by averaging all columns matching `"Biovolume"` and
#' then multiplying by the factor and cell count.
#'
#' @param df A BSA dataframe containing raw count and factor columns
#' @param unit_col Name of the column containing unit abundance counts
#'   (default: `'Unit Abundance'`)
#' @param cell_col Name of the column containing total cell counts
#'   (default: `'Total Number of Cells'`)
#' @param calc_cols Character vector specifying which density columns to
#'   calculate. One or more of `'Units'`, `'Cells'`, `'Biovolume'`
#'   (default: all three)
#'
#' @return Dataframe with one or more new density columns added:
#'   `Units_per_mL`, `Cells_per_mL`, `Biovolume_per_mL`
#'
#' @importFrom dplyr mutate across all_of select matches
#' @export
calc_data_bsa <- function(df,
                          unit_col = 'Unit Abundance',
                          cell_col = 'Total Number of Cells',
                          calc_cols = c('Units', 'Cells', 'Biovolume')) {

  calc_cols <- match.arg(calc_cols, choices = c('Units', 'Cells', 'Biovolume'), several.ok = TRUE)
  calculated <- character()

  # find factor column
  factor_col <- names(df)[grepl('factor', names(df), ignore.case = TRUE)]
  if (length(factor_col) == 0) {
    stop('Factor column not found')
  }

  # convert relevant cols to numeric
  df <- df %>%
    mutate(across(all_of(c(unit_col, cell_col)), as.numeric),
                  across(matches('Biovolume', ignore.case = TRUE), as.numeric),
                  Factor = as.numeric(.data[[factor_col]]))

  if ('Units' %in% calc_cols) {
    df <- df %>%
      mutate(Units_per_mL = round(.data[[unit_col]] * .data[[factor_col]], 2))
    calculated <- c(calculated, 'Units_per_mL')
  }

  if ('Cells' %in% calc_cols) {
    df <- df %>%
      mutate(Cells_per_mL = round(.data[[cell_col]] * .data[[factor_col]], 2))
    calculated <- c(calculated, 'Cells_per_mL')
  }

  # biovol is calc'd by averaging biovol columns
  if ('Biovolume' %in% calc_cols) {
    df <- df %>%
      mutate(Biovolume_per_mL = round(
        rowMeans(select(., matches('Biovolume', ignore.case = TRUE)), na.rm = TRUE) * .data[[factor_col]] * .data[[cell_col]],
        2
      ))
    calculated <- c(calculated, 'Biovolume_per_mL')
  }

  message('Calculated: ', paste(calculated, collapse = ', '))

  return(df)
}

#' @title Clean BSA PhytoForm codes
#'
#' @description
#' Standardizes the `PhytoForm` column by converting values to lowercase and
#' replacing non-informative entries (`"n/p"`, `"na"`, `"#n/a"`, `"0"`, `"1"`)
#' with `NA`. Prints a message listing the unique codes remaining after cleaning.
#'
#' @param df A dataframe containing a `PhytoForm` column
#'
#' @return Dataframe with a cleaned `PhytoForm` column
#'
#' @export
clean_phytoform_bsa <- function(df) {
  # lowercase letters
  df$PhytoForm <- tolower(df$PhytoForm)

  # replace NAs
  df$PhytoForm <- replace(df$PhytoForm, df$PhytoForm %in% c('n/p', 'na', 'NA', '#n/a', '0', '1'), NA_character_)

  unique_codes <- unique(df$PhytoForm)
  message('PhytoForm codes: ', paste(sort(unique_codes[!is.na(unique_codes)]), collapse = ', '))

  return(df)
}

# Main Function -----------------------------------------------------------

# Full BSA standardization pipeline: removes blank rows, renames columns,
# calculates density values, cleans PhytoForm, joins sampling method from
# metadata, and selects final output columns
#' @importFrom dplyr rename select all_of
#' @noRd
standardize_cols_bsa <- function(df, meta_df){
  df <- remove_rows_bsa(df)

  df <- rename_cols_bsa(df)

  df <- calc_units_bsa(df)

  df <- clean_phytoform_bsa(df)

  # pull collection type
  df <- from_meta(df, meta_df, 'Sampling Method')

  df <- df %>% rename('SamplingMethod' = `Sampling Method`)

  # select and reorder relevant columns
  df <- df %>%
    select(all_of(c('Lab', 'Date', 'Time', 'Station', 'Taxon', 'Genus', 'Species', 'Units_per_mL', 'Cells_per_mL', 'Biovolume_per_mL', 'GALD', 'PhytoForm', 'Comments')))
}


# EA Calcs ----------------------------------------------------------------

#' @title Calculate per-mL density values from EA data
#'
#' @description
#' Computes standardized density columns (`Units_per_mL`, `Cells_per_mL`,
#' and/or `Biovolume_per_mL`) from EA-formatted data.
#' `Units_per_mL` is taken directly from the unit column.
#' Cells are derived by multiplying units by cells-per-unit.
#' Biovolume requires a factor column and is calculated by averaging all columns
#' matching `"Biovolume"` and multiplying by factor and cell count.
#'
#' @param df An EA dataframe containing raw count columns
#' @param unit_col Name of the column containing organisms per mL
#'   (default: `'Organisms per mL'`)
#' @param cell_col Name of the column containing cells per unit
#'   (default: `'Cells per Unit'`)
#' @param calc_cols Character vector specifying which density columns to
#'   calculate. One or more of `'Units'`, `'Cells'`, `'Biovolume'`
#'   (default: all three)
#'
#' @return Dataframe with one or more new density columns added:
#'   `Units_per_mL`, `Cells_per_mL`, `Biovolume_per_mL`
#'
#' @importFrom dplyr mutate across any_of select matches
#' @export
calc_data_ea <- function(df,
                         unit_col = 'Organisms per mL',
                         cell_col = 'Cells per Unit',
                         calc_cols = c('Units', 'Cells', 'Biovolume')) {

  calc_cols <- match.arg(calc_cols, choices = c('Units', 'Cells', 'Biovolume'), several.ok = TRUE)
  calculated <- character()

  # only require factor column if Biovolume is calculated
  if ('Biovolume' %in% calc_cols) {
    factor_col <- names(df)[grepl('factor', names(df), ignore.case = TRUE)]
    if (length(factor_col) == 0) stop('Factor column not found')
  } else {
    factor_col <- NULL
  }

  # convert only existing numeric columns
  df <- df %>%
    mutate(across(any_of(c(unit_col, cell_col)), as.numeric),
           across(matches('Biovolume', ignore.case = TRUE), as.numeric))

  if (!is.null(factor_col))
    df <- df %>% mutate(Factor = as.numeric(.data[[factor_col]]))

  if ('Units' %in% calc_cols) {
    df <- df %>%
      mutate(Units_per_mL = round(.data[[unit_col]], 2))
    calculated <- c(calculated, 'Units_per_mL')
  }

  if ('Cells' %in% calc_cols && cell_col %in% names(df)) {
    df <- df %>%
      mutate(Cells_per_mL = round(.data[[unit_col]] * .data[[cell_col]], 2))
    calculated <- c(calculated, 'Cells_per_mL')
  }

  if ('Biovolume' %in% calc_cols && !is.null(factor_col)) {
    df <- df %>%
      mutate(Biovolume_per_mL = round(
        rowMeans(select(., matches('Biovolume', ignore.case = TRUE)), na.rm = TRUE) *
          .data[[factor_col]] * .data[[cell_col]],
        2
      ))
    calculated <- c(calculated, 'Biovolume_per_mL')
  }

  message('Calculated: ', paste(calculated, collapse = ', '))
  df
}
