
# Check Variables ---------------------------------------------------------

#' @title Check for non-distinct rows under different key definitions
#'
#' @description
#' Checks whether rows in a dataframe are distinct under one of three definitions:
#' the full row (`type = 'full'`), the combination of key + taxa + measurement columns
#' (`type = 'key_measure_cols'`), or the combination of key + taxa columns (`type = 'key_cols'`).
#'
#' If non-distinct rows are found, they are collected and either returned directly
#' (`return_df = TRUE`) or attached to the input dataframe via the `log` attribute.
#'
#' @param df input dataframe to check
#' @param return_df logical; if `TRUE`, returns the non-distinct rows as a dataframe
#'   instead of returning the original dataframe with a `log` attribute
#' @param coerce logical; if `TRUE`, coerces factor columns to character and rounds
#'   numeric columns to 6 decimals prior to distinctness checks
#' @param type character; definition of distinctness to check, one of
#'   `'full'`, `'key_measure_cols'`, or `'key_cols'`
#' @param key_cols character vector of key columns used to identify samples
#'   (default `c('Station','Date')`)
#' @param taxa_cols character vector of taxa identifier columns (default `c('Taxon')`)
#' @param measurement_cols character vector of measurement columns used when
#'   `type = 'key_measure_cols'` (default `c('Biovolume_per_mL','Cells_per_mL','Units_per_mL')`)
#'
#' @return
#' If no non-distinct rows are found, returns `invisible(df)` and stores `NULL` in the
#' relevant log entry.
#'
#' If non-distinct rows are found:
#' \itemize{
#'   \item when `return_df = TRUE`, returns a dataframe of non-distinct rows
#'   \item otherwise returns `df` with `attr(df, 'log')` updated
#' }
#'
#' The `log` attribute uses one of the following names depending on `type`:
#' \itemize{
#'   \item `nondistinct_allrows` for `type = 'full'`
#'   \item `nondistinct_keymeasurerows` for `type = 'key_measure_cols'`
#'   \item `nondistinct_keyrows` for `type = 'key_cols'`
#' }
#'
#' @importFrom dplyr mutate across where distinct group_by summarize filter select semi_join all_of
#' @export
check_distinct <- function(df, return_df = FALSE, coerce = TRUE,
                           type = c('full', 'key_measure_cols', 'key_cols'),
                           key_cols = c('Station','Date'),
                           taxa_cols = c('Taxon'),
                           measurement_cols = c('Biovolume_per_mL','Cells_per_mL','Units_per_mL')) {
  type <- match.arg(type)
  df_check <- df
  
  if (coerce) {
    df_check <- df_check %>%
      mutate(
        across(where(is.factor), as.character),
        across(where(is.numeric), ~ round(.x, 6))  
      )
  }
  
  if (type == 'full') {
    n_total <- nrow(df_check)
    n_unique <- nrow(distinct(df_check))
    n_non_distinct <- n_total - n_unique
    
    if (n_non_distinct > 0) {
      df_dupes <- df[duplicated(df_check) | duplicated(df_check, fromLast = TRUE), ]
      message('Number of non-distinct rows: ', nrow(df_dupes))
      
      if (isTRUE(return_df)) {
        return(df_dupes)
      } else {
        attr(df, 'log') <- list(nondistinct_allrows = df_dupes)
        return(df)
      }
    } else {
      message('All rows are unique.')
      attr(df, 'log') <- list(nondistinct_allrows = NULL)
      return(invisible(df))
    }
    
  } else if (type == 'key_measure_cols') {
    if (is.null(key_cols) || is.null(measurement_cols)) {
      stop('For type = "key_measure_cols", both key_cols and measurement_cols must be provided.')
    }
    
    # find non-distinct keys
    nondistinct_keys <- df_check %>%
      group_by(across(all_of(c(key_cols, taxa_cols, measurement_cols)))) %>%
      summarise(n = n(), .groups = 'drop') %>%
      filter(n > 1) %>%
      select(all_of(c(key_cols, taxa_cols, measurement_cols)))
    
    # extract all original rows that match these non-distinct keys
    df_dupes <- df %>%
      semi_join(nondistinct_keys, by = c(key_cols, taxa_cols, measurement_cols))
    
    n_inconsistent <- nrow(df_dupes)
    
    if (n_inconsistent > 0) {
      message('Number of inconsistent key+taxa+measurement rows: ', n_inconsistent)
      attr(df, 'log') <- list(nondistinct_keymeasurerows = df_dupes)
      if (isTRUE(return_df)) {
        return(df_dupes)
      } else {
        return(df)
      }
    } else {
      message('All key+taxa+measurement column combinations are unique.')
      attr(df, 'log') <- list(nondistinct_keymeasurerows = NULL)
      return(invisible(df))
    }
    
  } else if (type == 'key_cols') {
    # find non-distinct keys
    nondistinct_keys <- df_check %>%
      group_by(across(all_of(c(key_cols, taxa_cols)))) %>%
      summarise(n = n(), .groups = 'drop') %>%
      filter(n > 1) %>%
      select(all_of(c(key_cols, taxa_cols)))
    
    # extract all original rows that match these non-distinct keys
    df_dupes <- df %>%
      semi_join(nondistinct_keys, by = c(key_cols, taxa_cols))
    
    n_non_distinct <- nrow(df_dupes)
    
    if (n_non_distinct > 0) {
      message('Number of non-distinct key+taxa rows: ', n_non_distinct)
      attr(df, 'log') <- list(nondistinct_keyrows = df_dupes)
      if (isTRUE(return_df)) {
        return(df_dupes)
      } else {
        return(df)
      }
    } else {
      message('All key+taxa combinations are unique.')
      attr(df, 'log') <- list(nondistinct_keyrows = NULL)
      return(invisible(df))
    }
  }
}

#' @title Extract unstandardized comments from a dataframe
#'
#' @description
#' Identifies and extracts comment text fragments that do not match a predefined
#' set of known quality or debris phrases. Normalizes delimiters, removes known
#' phrases, and appends unmatched comment fragments to the dataframe's log
#' attribute under `unmatched_comments`.
#'
#' @param df a dataframe containing comment text
#' @param comment_col unquoted column name containing comments to parse
#' @param delimiter character string giving the delimiter used to separate comment fragments (default `' '`)
#'
#' @return the input dataframe with an updated `'log'` attribute containing a tibble of unmatched comment fragments (if any)
#'
#' @importFrom rlang ensym
#' @importFrom dplyr pull case_when
#' @importFrom stringr str_remove_all str_replace_all str_split str_trim
#' @importFrom purrr map discard
#' @importFrom tibble tibble
#'
#' @export
extract_unstandardized_comments <- function(df, comment_col, delimiter = ' ') {
  comment_col <- ensym(comment_col)
  
  # normalize common delimiters to regex equivalents
  if (!is.null(delimiter) && delimiter != '') {
    delimiter <- case_when(
      delimiter == '. '  ~ '\\.\\s*',
      delimiter == '.'   ~ '\\.',
      delimiter == ', '  ~ ',\\s*',
      delimiter == ','   ~ ',',
      delimiter == '; '  ~ ';\\s*',
      delimiter == ';'   ~ ';',
      delimiter == ' - ' ~ '\\s+-\\s*',
      TRUE               ~ delimiter
    )
  }
  
  # phrases captured by QualityCheck or Debris
  known_phrases <- c(
    'delete',
    'cross contamination',
    'did not reach',
    'cannot meet tally',
    'cannot meet natural unit',
    'CNMT>5', 'CNMT >5', 'CNMT > 5',
    'CMT>5',  'CMT >5',  'CMT > 5',
    'CMNT<5', 'CMNT <5', 'CMNT < 5',
    'CMT<5',  'CMT <5',  'CMT < 5',
    'CMT', 'CMNT','LessThan400cells',
    'degraded',
    'poor preservation', 'poorly preserved', 'PoorlyPreserved',
    'weak preservation', 'weakly preserved',
    'fungus',
    'fungal growth',
    'mycelial growth',
    'obscured',
    'girdle','girdle view',
    'cyst',
    'mucilaginous detritus',
    'ciliates','ciliates observed',
    'many broken diatoms', 'broken diatoms', 'BrokenDiatoms',
    'high detritus', 'high sediment',
    'moderate detritus', 'moderate sediment','moderat sediment',
    'low detritus', 'low sediment',
    'light detritus', 'light sediment',
    'heavy detritus', 'heavy sediment', 'high amount of debris',
    'high amounts of debris', 'lots of debris', 'High sedimnet', 'High sedimen',
    'fragment.', 'fragmented diatoms', 'diatom fragments', 'diatom fragment',
    'Good',
    'and'
  )
  
  # compile regex pattern
  known_pattern <- paste0(
    '(', 
    paste0(str_replace_all(known_phrases, '(\\W)', '\\\\\\1'), collapse = '|'), 
    ')'
  )
  
  raw_comments <- df %>%
    pull(!!comment_col) %>%
    unique() %>%
    discard(is.na)
  
  cleaned_comments <- str_remove_all(raw_comments, regex(known_pattern, ignore_case = TRUE))
  
  # isolate unmatched fragments
  if (!is.null(delimiter) && delimiter != '') {
    unmatched <- cleaned_comments %>%
      map(~ str_split(.x, delimiter)[[1]]) %>%
      unlist() %>%
      str_trim() %>%
      str_remove('\\.$') %>%
      str_remove('\\;$') %>%
      discard(~ .x == '' || is.na(.x))
  } else {
    unmatched <- cleaned_comments %>%
      str_trim() %>%
      discard(~ .x == '' || is.na(.x))
  }
  
  unmatched_df <- unique(unmatched) %>%
    tibble(Unmatched = .)
  
  if (nrow(unmatched_df) > 0) {
    message('Unique unstandardized comments: ', nrow(unmatched_df))
    attr(df, 'log')$unmatched_comments <- unmatched_df
  }
  
  return(df)
}

# Check Taxa --------------------------------------------------------------

#' @title Check for missing higher-level taxonomic classifications
#'
#' @description
#' Extracts the taxonomic classification columns from `Taxon` through
#' `Species` and identifies taxa that are missing one or more higher-level
#' fields (`Kingdom`, `Phylum`, `Class`, `AlgalGroup`,
#' `Genus`, `Species`). Prior to reporting, the function standardizes
#' `Taxon` by removing `cf.` and trimming variety suffixes
#' (anything matching `' var\\..*`), then de-duplicates and sorts the result.
#'
#' If any taxa are missing classifications, a warning is issued listing the unique
#' affected taxa. Otherwise, a message is printed indicating classifications are present.
#'
#' @param df input dataframe containing columns `Taxon` through `Species`
#'
#' @return A dataframe containing the rows for taxa that are missing
#' one or more higher-level classification fields.
#' Returns an empty dataframe if none are missing.
#'
#' @importFrom dplyr select mutate arrange
#' @importFrom glue glue
#' @export
check_higher_taxa <- function(df){
  df_error_check <- df %>% select(c('Taxon':'Species'))
  
  check <- df_error_check %>% subset(is.na(AlgalGroup) | is.na(Class) | is.na(Phylum) | is.na(Kingdom) | is.na(Genus) | is.na(Species))
  check <- check %>% mutate(Taxon = gsub('cf\\. ', '', Taxon),
                            Taxon = gsub(' var\\..*', '', Taxon))
  check <- check[!duplicated(check),]
  check <- check %>% arrange(Taxon)  
  
  if (nrow(check) > 0){
    warning(glue('Taxon missing higher level classifications:\n{toString(unique(check$Taxon)\n)}\nEither update official list or fix name(s)'))
  } else {
    message('All higher level classifications added.')
  }
  
  return(check)
}

#' @title Check synonym reference table for multi-generational chains and flag missing synonym entries
#'
#' @description
#' Reads a synonym/classification reference table and performs two checks:
#' \itemize{
#'   \item Identifies multi-generational synonym chains in the reference table,
#'   where a taxon that is itself updated (has `CurrentTaxon != 'None'`) also
#'   appears as a `CurrentTaxon` for some other record
#'   \item Identifies taxa in the input dataframe that are missing synonym coverage,
#'   defined here as rows where `Taxon` is `NA` or equals `'Unknown'`
#' }
#'
#' The function prints a message that includes any multi-generational warning text
#' and, if missing synonym taxa are present, lists the unique affected taxa.
#'
#' @param df input dataframe containing columns `Taxon` through `Species`
#'
#' @return A dataframe combining:
#' \itemize{
#'   \item a summary table of multi-generational taxa from the synonym reference
#'   \item the `Taxon` through `Species` rows from `df` where
#'   `Taxon` is missing (`NA`) or equals `Unknown`
#' }
#' The returned dataframe includes a `Type` column with values `multigen`
#' and `missing syn` to distinguish record types.
#'
#' @importFrom dplyr select mutate arrange full_join
#' @importFrom glue glue
#' @export
check_synonyms <- function(df){
  df_syn <- read_quiet_csv('admin/global_data/phyto_classifications.csv') %>%
    select(c('Kingdom':'AlgalGroup','Taxon','CurrentTaxon'))
  
  changed_taxon <- unique(df_syn$Taxon[df_syn$CurrentTaxon != 'None'])
  
  multigen <- changed_taxon[changed_taxon %in% unique(df_syn$CurrentTaxon)]
  
  df_output <- data.frame(multigen) %>%
    mutate(Type = 'multigen')
  
  if (length(multigen) > 0){
    warning_one <- glue('Warning: multi-generational taxon names in synonym dataframe: {toString(multigen)}.')
  } else {
    warning_one <- ''
  }
  
  df_error_check <- df %>% select(c('Taxon':'Species'))
  
  syn_check <- df_error_check %>% subset(is.na(Taxon) | Taxon == 'Unknown')
  syn_check <- syn_check[!duplicated(syn_check),]
  syn_check <- syn_check %>% arrange(Taxon)
  
  if (nrow(syn_check) > 0){
    message(glue('{warning_one} \n Warning: Taxon missing synonym data: {toString(unique(syn_check$Taxon)\n)}\n Either update official list or fix name(s)'))
  } else {
    message(glue('{warning_one} \n All synonyms added.'))
  }
  
  syn_check <- syn_check %>%
    mutate(Type = 'missing syn')
  
  df_output <- full_join(df_output, syn_check, by = 'Type')
  
  return(df_output)
}


# Check data --------------------------------------------------------------

#' @title Check for "No organisms observed" records
#'
#' @description
#' Filters the input dataframe for rows where `Taxon` equals
#' `'No organisms observed'` and prints a summary message listing the
#' associated `Station` and `Date` values. If no such rows are present,
#' prints `'No data observed for: None'`.
#'
#' @param df input dataframe containing at least `Taxon`, `Station`, and `Date`
#'
#' @return Invisibly returns a dataframe of the "no organisms observed" rows.
#' If none are found, returns an empty dataframe.
#'
#' @importFrom dplyr filter
#' @export
check_nodata <- function(df) {
  df_nodat <- df %>% filter(Taxon == 'No organisms observed')
  if (nrow(df_nodat) == 0) {
    print('No data observed for: None')
  } else {
    message <- paste0(
      'No data observed for: ',
      paste(df_nodat$Station, df_nodat$Date, sep = ' ', collapse = '; ')
    )
    print(message)
  }
}

#' @title Check that Cells_per_mL >= Units_per_mL
#'
#' @description
#' Verifies that `Cells_per_mL` is greater than or equal to `Units_per_mL` for
#' every row. Rows where cells are less than units are flagged and attached to
#' the dataframe's `log` attribute under `cell_calc_issue`. If either column is
#' absent, a message is printed and the dataframe is returned unchanged.
#'
#' @param df A dataframe containing `Units_per_mL` and `Cells_per_mL` columns
#'
#' @return The input dataframe, with a `log` attribute entry `cell_calc_issue`
#'   added if any rows fail the check
#'
#' @importFrom dplyr filter
#' @export
check_units_cells <- function(df) {
  if ('Units_per_mL' %in% colnames(df) && 'Cells_per_mL' %in% colnames(df)) {
    
    # Check that Cells_per_mL >= Units_per_mL for every row pair
    count_issues <- df %>%
      filter(Cells_per_mL < Units_per_mL)
    
    if (nrow(count_issues) == 0) {
      message('(Correct) Cells_per_mL is greater than or equal to Units_per_mL for all rows.')
    } else {
      log_df <- count_issues
      message('(Warning) Found ', nrow(count_issues), ' row where Cells_per_mL is less than Units_per_mL.')
      
      # # Add the QC code 'UnitsExceedCells' to the QualityCheck column
      # df <- df %>%
      #   mutate(QualityCheck = case_when(
      #     row_number() %in% rownames(count_issues) & QualityCheck != 'NoCode' ~ paste(QualityCheck, 'UnitsExceedCells'),
      #     row_number() %in% rownames(count_issues) & QualityCheck == 'NoCode' ~ 'UnitsExceedCells',
      #     TRUE ~ QualityCheck
      #   ))
      
      attr(df, 'log') <- list(cell_calc_issue = log_df)
    }
  } else {
    message('One or both of the required columns ("Units_per_mL", "Cells_per_mL") are missing.')
  }
  
  return(df)
}

#' @title Check for missing values across columns
#'
#' @description
#' Scans a dataframe for columns containing any `NA` values (partial NAs) and
#' columns where all values are `NA` (fully missing). Specified columns can be
#' excluded from the partial NA check. Results are printed as messages and
#' attached to the dataframe's `log` attribute under `na_check`.
#'
#' @param df A dataframe to check for missing values
#' @param exclude_cols Character vector of column names to exclude from the
#'   partial NA check (default: `'OrigTaxon'`)
#'
#' @return The input dataframe with a `log` attribute entry `na_check`
#'   containing a tibble of column names with partial NAs
#'
#' @importFrom tibble tibble
#' @importFrom dplyr distinct
#' @export
check_nas <- function(df, exclude_cols = 'OrigTaxon') {
  # Identify columns to check for partial NAs
  cols_to_check_partial <- if (!is.null(exclude_cols)) {
    setdiff(names(df), exclude_cols)
  } else {
    names(df)
  }
  
  # Check for partial NAs
  df_check_partial <- df[, cols_to_check_partial, drop = FALSE]
  na_cols_partial <- names(which(colSums(is.na(df_check_partial)) > 0))
  
  # Check for fully missing columns (including excluded columns)
  na_cols_full <- names(which(colSums(is.na(df)) == nrow(df)))
  
  # Create log as a tibble
  log_df <- tibble(MissingInColumn = na_cols_partial) %>% distinct()
  
  # Print messages
  if (length(na_cols_partial) > 0) {
    message('Missing values found in ', length(na_cols_partial), ' column(s): ', paste(na_cols_partial, collapse = ', '))
  } else {
    message('No missing values found (excluding ', paste(exclude_cols, collapse = ', '), ').')
  }
  
  if (length(na_cols_full) > 0) {
    message('All values are missing in ', length(na_cols_full), ' column(s): ', paste(na_cols_full, collapse = ', '))
  } else {
    message('No column has all missing values.')
  }
  
  # Attach log to df
  attr(df, 'log') <- list(na_check = log_df)
  
  return(df)
}

#' @title Print unique values in a column
#'
#' @description
#' Prints a message listing all unique values found in the specified column.
#'
#' @param df A dataframe
#' @param col Unquoted column name to inspect
#'
#' @return Called for its side effect (message); returns `NULL` invisibly
#'
#' @importFrom rlang ensym as_string
#' @importFrom dplyr pull
#' @export
unique_check <- function(df, col) {
  col_sym <- ensym(col)
  unique_vals <- unique(df %>% pull(!!col_sym))
  message('Unique ', as_string(col_sym), 's: ', paste(unique_vals, collapse = ', '))
}

#' @title Plot number of stations sampled per month
#'
#' @description
#' Summarizes the number of distinct stations sampled in each month and year,
#' then displays the result as a heatmap tile plot. Months with no sampling
#' are shown as white tiles.
#'
#' @param df A dataframe containing `Date` and `Station` columns
#'
#' @return A `ggplot` object
#'
#' @importFrom dplyr mutate distinct count left_join
#' @importFrom lubridate floor_date year month
#' @importFrom tidyr expand_grid
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_distiller scale_y_continuous theme_minimal labs theme element_text element_blank
#' @export
check_station_count <- function(df) {
  df <- df %>%
    mutate(
      MonthYear = as.Date(floor_date(Date, 'month')),
      Station = as.character(Station)
    )
  df_station_count <- df %>%
    distinct(Station, MonthYear) %>%
    count(MonthYear, name = 'NumStations') %>%
    mutate(
      Year = year(MonthYear),
      Month = month(MonthYear, label = TRUE, abbr = TRUE)
    )
  
  # Create complete grid for all year-month combinations
  all_years <- min(df_station_count$Year):max(df_station_count$Year)
  all_months <- month(1:12, label = TRUE, abbr = TRUE)  # This creates ordered factor like df_station_count$Month
  complete_grid <- expand_grid(Year = all_years, Month = all_months)
  
  df_complete <- complete_grid %>%
    left_join(df_station_count, by = c('Year', 'Month'))
  
  ggplot(df_complete, aes(x = Month, y = Year, fill = NumStations)) +
    geom_tile(color = 'black', linewidth = 0.5) +
    scale_fill_distiller(palette = 'Blues', direction = 1, name = '# Stations', na.value = 'white') +
    scale_y_continuous(breaks = all_years) +
    theme_minimal() +
    labs(
      x = 'Month',
      y = 'Year',
      title = 'Number of Stations Sampled Per Month'
    ) +
    theme(
      axis.text.x = element_text(angle = 45),
      panel.grid = element_blank()
    )
}

#' @title Plot station sampling frequency per year
#'
#' @description
#' Summarizes how many times each station was sampled per year, then displays
#' the result as a heatmap tile plot. Station-year combinations with no
#' sampling are shown as white tiles.
#'
#' @param df A dataframe containing `Date` and `Station` columns
#'
#' @return A `ggplot` object
#'
#' @importFrom dplyr mutate distinct count left_join
#' @importFrom lubridate floor_date year
#' @importFrom tidyr expand_grid
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_distiller scale_x_continuous theme_minimal labs theme element_text element_blank
#' @export
check_sampling_freq <- function(df) {
  df <- df %>%
    mutate(
      MonthYear = as.Date(floor_date(Date, 'month')),
      Station = as.character(Station)
    )
  
  df_station_frequency <- df %>%
    mutate(Year = year(MonthYear)) %>%
    distinct(Station, MonthYear, Year) %>%
    count(Station, Year, name = 'TimesPerYear')
  
  # Create complete grid for all station-year combinations
  all_years <- min(df_station_frequency$Year):max(df_station_frequency$Year)
  all_stations <- unique(df$Station)
  complete_grid <- expand_grid(Year = all_years, Station = all_stations)
  
  df_complete <- complete_grid %>%
    left_join(df_station_frequency, by = c('Year', 'Station'))
  
  ggplot(df_complete, aes(x = Year, y = Station, fill = TimesPerYear)) +
    geom_tile(color = 'black', linewidth = 0.2) +
    scale_fill_distiller(palette = 'Blues', direction = 1, name = 'Frequency', na.value = 'white') +
    scale_x_continuous(breaks = all_years) +
    theme_minimal() +
    labs(
      x = 'Year',
      y = 'Station',
      title = 'Station Sampling Frequency Per Year'
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
}

# Check Introduced NAs ----------------------------------------------------

#' @title Identify values that become NA when coerced to numeric
#'
#' @description
#' Attempts to coerce the specified column to numeric and reports any values
#' that become `NA` as a result (i.e. non-numeric strings that would be silently
#' lost in a conversion). Prints the problematic values and the row indices
#' where they occur (up to 10).
#'
#' @param df A dataframe containing the column to inspect
#' @param col_name Name of the column to check (as a character string)
#'
#' @return Called for its side effect (printed output); returns `NULL` invisibly
#'
#' @export
check_non_numeric <- function(df, col_name) {
  if (col_name %in% names(df)) {
    original <- df[[col_name]]
    converted <- as.numeric(as.character(original))
    
    # Find which values became NA
    na_indices <- which(is.na(converted) & !is.na(original))
    
    if (length(na_indices) > 0) {
      cat('Column:', col_name, '\n')
      cat('Problematic values:\n')
      print(unique(original[na_indices]))
      cat('At rows:', na_indices[1:min(10, length(na_indices))], '\n\n')
    }
  }
}