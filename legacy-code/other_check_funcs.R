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