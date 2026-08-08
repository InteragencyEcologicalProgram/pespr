#' @title Flag statistical outliers in a measurement column
#'
#' @description
#' Identifies outliers in log-transformed measurement data using robust z-scores 
#' (median and MAD), flags them in the \code{QualityCheck} column, and optionally 
#' generates a scatter plot showing flagged points. The function also appends a 
#' dataframe of flagged rows to the \code{log} attribute.
#'
#' @param df input dataframe containing the specified measurement column and 
#'   \code{QualityCheck}, \code{Notes}, \code{Date}, \code{Station}, and \code{Taxon} columns
#' @param col unquoted column name of the measurement variable to evaluate
#'   (must be one of \code{Cells_per_mL}, \code{Units_per_mL}, or
#'   \code{Biovolume_per_mL})
#' @param station_col unquoted column name for the station grouping variable (default = \code{Station})
#' @param cutoff numeric value for the absolute robust z-score threshold
#'   used to identify outliers (default = 3)
#' @param add_flag logical; if \code{TRUE}, updates \code{QualityCheck} for flagged rows (default = TRUE)
#' @param show_plot logical; if \code{TRUE}, prints a scatter plot with
#'   flagged outliers in red (default = TRUE)
#'   
#' @return the input dataframe with updated \code{QualityCheck} values and 
#' an updated \code{log} attribute containing a dataframe of flagged rows
#'   
#' @importFrom rlang enquo as_name
#' @importFrom dplyr mutate filter pull select case_when na_if
#' @importFrom ggplot2 ggplot aes geom_point scale_fill_manual scale_color_manual scale_y_log10 labs theme_minimal
#' @importFrom scales alpha
#' @importFrom stringr str_detect
#' @importFrom stats median mad
#' @importFrom lubridate year
#' @export
flag_outliers <- function(df, col, station_col = Station, cutoff = 3, add_flag = TRUE, show_plot = TRUE) {
  col <- enquo(col)
  col_name <- as_name(col)
  station <- enquo(station_col)
  station_name <- as_name(station)
  
  # determine code label based on column name (case-insensitive, flexible)
  switch_label <- case_when(
    str_detect(tolower(col_name), 'cells') ~ 'Cells',
    str_detect(tolower(col_name), 'units') ~ 'Units',
    str_detect(tolower(col_name), 'biovol') ~ 'Biovol',
    TRUE ~ NA_character_
  )
  if (is.na(switch_label)) {
    stop('colname must contain one of: "cells", "units", or "biovolume"')
  }
  
  # compute robust z-scores (upper outliers only)
  df_flagged <- df %>%
    mutate(
      log_val = if_else(!!col > 0, log10(!!col), NA_real_),
      z_robust = 0.6745 * (log_val - median(log_val, na.rm = TRUE)) /
        mad(log_val, na.rm = TRUE),
      outlier = if_else(!is.na(z_robust) & z_robust > cutoff, TRUE, FALSE)
    )
  
  n_out <- sum(df_flagged$outlier, na.rm = TRUE)
  message(n_out, ' outlier(s) flagged in ', col_name)
  
  # compute percentage above
  pct_above <- mean(df_flagged$z_robust > cutoff, na.rm = TRUE)
  if (pct_above > 0) {
    message(sprintf('%.3f%% of data above %.2f MAD threshold', pct_above * 100, cutoff))
  } else {
    message(sprintf('No data above %.2f MAD threshold', cutoff))
  }
  
  df_out <- df
  
  if(add_flag){
    # update QualityCheck column
    qc_update <- df_flagged %>%
      mutate(
        QualityCheck_clean = str_remove_all(QualityCheck, paste0('\\s*;?\\s*Outlier', switch_label)) %>%
          str_squish() %>%
          na_if(''),
        QualityCheck_clean = case_when(
          is.na(QualityCheck_clean) ~ 'NoCode',
          TRUE ~ QualityCheck_clean
        ),
        QualityCheck_new = case_when(
          outlier & QualityCheck_clean == 'NoCode' ~ paste0('Outlier', switch_label),
          outlier & QualityCheck_clean != 'NoCode' ~ paste(QualityCheck_clean, paste0('Outlier', switch_label), sep = '; '),
          TRUE ~ QualityCheck_clean
        ) %>%
          gsub('^;\\s*|NA;\\s*', '', .)
      ) %>%
      pull(QualityCheck_new)
    
    # overwrite QualityCheck in original df
    df_out$QualityCheck <- qc_update
  } 
  
  # optional plot
  if (show_plot) {
    p <- df_flagged %>%
      filter(!is.na(!!col)) %>%
      mutate(x_index = cumsum(!is.na(!!col))) %>%
      ggplot(aes(x = x_index, y = !!col)) +
      geom_point(aes(fill = outlier), color = 'black', shape = 21, size = 2.5) +
      scale_fill_manual(values = c('FALSE' = 'white', 'TRUE' = 'red')) +
      labs(
        y = col_name,
        x = NULL,
        title = paste('Flagged outliers in', col_name)
      ) +
      theme_minimal()
    print(p)
  }
  
  # add log of flagged rows
  existing_log <- attr(df, 'log')
  outlier_rows <- df_out %>%
    filter(str_detect(QualityCheck, paste0('Outlier', switch_label))) %>%
    select(Date, !!station, Taxon, !!col, QualityCheck, any_of('Notes'))
  
  attr(df_out, 'log') <- c(
    existing_log,
    setNames(list(outlier_rows), paste0('outlier_', tolower(switch_label)))
  )
  
  return(df_out)
}