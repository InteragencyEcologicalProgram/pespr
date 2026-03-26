
# Plot NMDS ---------------------------------------------------------------

#' @title Build a community matrix from long-format data
#'
#' @description
#' Aggregates a long-format dataframe into a wide community matrix suitable
#' for multivariate analysis. Optionally applies an effort correction by
#' dividing totals by the number of unique sampling events in `avg_var`.
#' Returns both the full combined dataframe and the raw community matrix.
#'
#' @param df A long-format dataframe
#' @param group_var Character vector of column names defining each sample event
#'   (each unique combination becomes one row in the matrix)
#' @param measure_var Name of the numeric column to aggregate (e.g. `'Cells_per_mL'`)
#' @param taxa_var Name of the taxon column to pivot to wide format
#'   (default: `'Taxon'`)
#' @param factor_var Optional character vector of additional metadata columns
#'   to retain (e.g. grouping variables for coloring plots)
#' @param avg_var Optional character vector of columns to use for effort
#'   correction; totals are divided by the number of distinct `avg_var` events
#'   within each `group_var` combination
#' @param compute_dist Logical; if `TRUE`, computes a dissimilarity matrix
#'   using `vegdist()` and attaches it to the output as `m_dist`
#'   (default: `FALSE`)
#' @param distance Distance metric to pass to `vegdist()` (e.g. `'bray'`);
#'   only used when `compute_dist = TRUE`
#'
#' @return A named list with elements:
#' \itemize{
#'   \item `m_all` — combined metadata and community matrix as a dataframe
#'   \item `m_com` — numeric community matrix
#'   \item `m_meta` — metadata columns only
#'   \item `meta_vars` — list of the variable name arguments
#'   \item `m_dist` — dissimilarity matrix (only present if `compute_dist = TRUE`)
#' }
#'
#' @importFrom dplyr distinct across all_of count group_by summarise left_join mutate select rename bind_cols
#' @importFrom tidyr pivot_wider
#' @export
create_mcom <- function(df, group_var, measure_var, taxa_var = 'Taxon',
                        factor_var = NULL, avg_var = NULL, compute_dist = FALSE, distance = NULL) {
  set.seed(42)

  message(
    '\n--- Params ---\n',
    'Measurement Variable: ', measure_var, '\n',
    'Multivariate Variable: ', taxa_var, '\n',
    'Event Grouping Variable: ', paste(group_var, collapse = ', '), ' (defines each point)\n',
    'Event Averaging Variable(s): ',
    if (!is.null(avg_var)) paste(avg_var, collapse = ', ') else 'none (grouping variable is at the sampling event level)', '\n',
    '\n---------------\n'
  )

  # Pre-compute grouping columns once
  group_cols <- c(group_var, factor_var, taxa_var)
  meta_vars <- c(group_var, factor_var)

  # Compute effort once if needed
  effort_divisor <- NULL
  if (!is.null(avg_var)) {
    effort_divisor <- df %>%
      distinct(across(all_of(c(group_var, avg_var)))) %>%
      count(across(all_of(group_var)), name = 'n_effort')
  }

  # aggregate
  df_sum <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(total = sum(.data[[measure_var]], na.rm = TRUE), .groups = 'drop')

  # Apply effort correction if needed
  if (!is.null(effort_divisor)) {
    df_sum <- df_sum %>%
      left_join(effort_divisor, by = group_var) %>%
      mutate(value = total / n_effort) %>%
      select(-total, -n_effort)
  } else {
    df_sum <- df_sum %>%
      rename(value = total)
  }

  # Pivot and prepare matrix (combined operation)
  df_cells <- df_sum %>%
    pivot_wider(
      names_from = all_of(taxa_var),
      values_from = value,
      values_fill = 0
    )

  # Extract metadata and community matrix in one step
  meta <- df_cells %>% select(all_of(meta_vars))
  m_com <- df_cells %>%
    select(-all_of(meta_vars)) %>%
    as.matrix()

  result <- bind_cols(
    meta,
    m_com
  )

  result_list <- list(
    m_all = result,
    m_com = m_com,
    m_meta = meta,
    meta_vars = list(
      measure_var = measure_var,
      taxa_var = taxa_var,
      group_var = group_var,
      factor_var = factor_var,
      avg_var = avg_var
    )
  )

  if (compute_dist) {
    dist_obj <- vegdist(m_com, method = distance)
    result_list$m_dist <- dist_obj
  }

  return(result_list)
}

#' @title Run NMDS on long-format community data
#'
#' @description
#' Aggregates long-format data into a community matrix and runs non-metric
#' multidimensional scaling (NMDS) via `vegan::metaMDS()`. Optionally applies
#' effort correction if `avg_var` is provided. Returns NMDS site scores joined
#' with sample metadata.
#'
#' @param df A long-format dataframe
#' @param group_var Character vector of column names defining each sample event
#'   (each unique combination becomes one point in the ordination)
#' @param measure_var Name of the numeric column to aggregate (e.g. `'Cells_per_mL'`)
#' @param taxa_var Name of the taxon column (default: `'Taxon'`)
#' @param factor_var Optional character vector of metadata columns to retain
#'   alongside NMDS scores (e.g. for coloring plots)
#' @param avg_var Optional character vector of columns for effort correction
#' @param distance Distance metric to use (default: `'bray'`)
#' @param show_legend Logical; currently unused, reserved for future use
#' @param color_palette Currently unused, reserved for future use
#' @param try Number of random starts for NMDS (default: `20`)
#' @param trymax Maximum number of random starts (default: `20`)
#'
#' @return A named list with elements:
#' \itemize{
#'   \item `df_nmds` — dataframe of NMDS1/NMDS2 site scores joined with metadata
#'   \item `m_com` — numeric community matrix used in the ordination
#'   \item `m_meta` — metadata columns only
#'   \item `meta_vars` — list of the variable name arguments
#' }
#'
#' @importFrom dplyr distinct across all_of count group_by summarise left_join mutate select rename bind_cols
#' @importFrom tidyr pivot_wider
#' @export
calc_nmds <- function(df, group_var, measure_var, taxa_var = 'Taxon',
                      factor_var = NULL, avg_var = NULL,
                      distance = 'bray',
                      show_legend = TRUE, color_palette = NULL,
                      try = 20, trymax = 20) {
  set.seed(42)

  message(
    '\n--- NMDS Setup ---\n',
    'Measurement Variable: ', measure_var, '\n',
    'Multivariate Variable: ', taxa_var, '\n',
    'Event Grouping Variable: ', group_var, ' (defines each point)\n',
    'Event Averaging Variable(s): ',
    if (!is.null(avg_var)) paste(avg_var, collapse = ', ') else 'none (grouping variable is at the sampling event level)', '\n',
    'Display/Coloring Variable(s): ',
    if (!is.null(factor_var)) paste(factor_var, collapse = ', ') else 'none',
    '\nDistance Metric: ', distance,
    '\n------------------\n'
  )

  # Pre-compute grouping columns once
  group_cols <- c(group_var, factor_var, taxa_var)
  meta_vars <- c(group_var, factor_var)

  # Compute effort once if needed
  effort_divisor <- NULL
  if (!is.null(avg_var)) {
    effort_divisor <- df %>%
      distinct(across(all_of(c(group_var, avg_var)))) %>%
      count(across(all_of(group_var)), name = 'n_effort')
  }

  # Single aggregation step
  df_sum <- df %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(total = sum(.data[[measure_var]], na.rm = TRUE), .groups = 'drop')

  # Apply effort correction if needed
  if (!is.null(effort_divisor)) {
    df_sum <- df_sum %>%
      left_join(effort_divisor, by = group_var) %>%
      mutate(value = total / n_effort) %>%
      select(-total, -n_effort)
  } else {
    df_sum <- df_sum %>%
      rename(value = total)
  }

  # Pivot and prepare matrix (combined operation)
  df_cells <- df_sum %>%
    pivot_wider(
      names_from = all_of(taxa_var),
      values_from = value,
      values_fill = 0
    )

  # Extract metadata and community matrix in one step
  meta <- df_cells %>% select(all_of(meta_vars))
  m_com <- df_cells %>%
    select(-all_of(meta_vars)) %>%
    as.matrix()

  # Run NMDS
  nmds <- suppressMessages(
    metaMDS(m_com, distance = distance, try = try, trymax = trymax)
  )

  # Combine results
  result <- bind_cols(
    as.data.frame(scores(nmds)$sites),
    meta
  )

  return(list(
    df_nmds = result,
    # raw_nmds = nmds,
    m_com = m_com,
    m_meta = meta,
    meta_vars = list(
      measure_var = measure_var,
      taxa_var = taxa_var,
      group_var = group_var,
      factor_var = factor_var,
      avg_var = avg_var
    )
  ))
}

#' @title Plot NMDS ordination results
#'
#' @description
#' Produces a scatter plot of NMDS1 vs NMDS2 site scores, colored by a
#' specified fill variable. Points are drawn as filled circles with a black
#' outline. Optionally applies a custom color palette and suppresses the legend.
#'
#' @param df_nmds A dataframe of NMDS site scores as returned by `calc_nmds()`
#' @param meta_vars List of variable name metadata as returned by `calc_nmds()`
#'   (used to auto-generate the plot title if `title` is `NULL`)
#' @param fill_var Name of the column to use for fill color (as a character string)
#' @param show_legend Logical; if `FALSE`, the legend is hidden (default: `TRUE`)
#' @param color_palette Optional named character vector of colors to pass to
#'   `scale_fill_manual()`
#' @param title Optional plot title; auto-generated from `meta_vars` if `NULL`
#' @param alpha Point transparency (default: `0.7`)
#' @param flip_order Logical; if `TRUE`, points are drawn in descending order
#'   of `fill_var` so lower values appear on top (default: `FALSE`)
#'
#' @return A `ggplot` object
#'
#' @importFrom rlang sym
#' @importFrom dplyr arrange desc
#' @importFrom ggplot2 ggplot aes geom_point theme_bw labs scale_fill_manual theme
#' @export
plot_nmds <- function(df_nmds, meta_vars, fill_var, show_legend = TRUE, color_palette = NULL, title = NULL, alpha = 0.7, flip_order = FALSE) {
  fill_sym <- rlang::sym(fill_var)

  if (is.null(title)) {
    title <- paste0(meta_vars$measure_var,
                    ' (per ', meta_vars$taxa_var,
                    ') per ', meta_vars$group_var)
  }

  if (flip_order) {
    df_nmds <- df_nmds %>%
      arrange(desc(!!fill_sym))
  } else {
    df_nmds <- df_nmds %>%
      arrange(!!fill_sym)
  }

  plt <- ggplot(df_nmds, aes(x = NMDS1, y = NMDS2, fill = !!fill_sym)) +
    geom_point(size = 4, shape = 21, alpha = alpha, color = '#000000') +
    theme_bw() +
    labs(title = title)

  if (!is.null(color_palette)) {
    plt <- plt + scale_fill_manual(values = color_palette)
  }

  if (!show_legend) {
    plt <- plt + theme(legend.position = 'none')
  }

  return(plt)
}

# Plot SIMPER Taxon -------------------------------------------------------

#' @title Plot top SIMPER taxa as boxplots
#'
#' @description
#' Extracts the top taxa from a SIMPER result object (filtered by significance
#' threshold `p`) and plots their distribution across `x_axis` groups as
#' side-by-side boxplots, one panel per taxon.
#'
#' @param df A dataframe containing the measurement and grouping columns
#' @param sim A SIMPER result object as returned by `vegan::simper()`
#' @param taxon_col Name of the taxon column (default: `'Taxon'`)
#' @param x_axis Name of the x-axis grouping column (default: `'Year'`)
#' @param y_axis Name of the y-axis measurement column
#' @param p Significance threshold for filtering SIMPER taxa (default: `0.05`)
#' @param taxa_num Maximum number of top taxa to display (default: `6`)
#'
#' @return A `ggplot` object
#'
#' @importFrom dplyr filter arrange pull
#' @importFrom ggplot2 ggplot aes geom_boxplot facet_wrap vars labs theme_minimal theme element_text
#' @importFrom tibble rownames_to_column
#' @export
plot_simper_bp <- function(df, sim, taxon_col = 'Taxon', x_axis = 'Year', y_axis, p = 0.05, taxa_num = 6) {
  taxon_sym <- sym(taxon_col)
  x_sym <- sym(x_axis)
  y_sym <- sym(y_axis)

  # extract top taxa from SIMPER results
  top_taxa <- sim[[1]] %>%
    as.data.frame() %>%
    rownames_to_column(var = taxon_col) %>%
    arrange(desc(average)) %>%
    filter(p <= !!p) %>%
    slice_head(n = taxa_num) %>%
    pull(!!taxon_sym)

  df %>%
    filter((!!taxon_sym) %in% top_taxa) %>%
    ggplot(aes(x = !!x_sym, y = !!y_sym)) +
    geom_boxplot(color = 'black') +
    facet_wrap(vars(!!taxon_sym), scales = 'free_y') +
    labs(
      title = paste('Top SIMPER', y_axis),
      x = x_axis,
      y = y_axis
    ) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8))
}

#' @title Plot top SIMPER taxa as aggregated line summaries
#'
#' @description
#' Extracts the top taxa from a SIMPER result object (filtered by significance
#' threshold `p`), aggregates their values over `x_axis` using `agg_fun`, and
#' plots the result as a line chart with one line per taxon.
#'
#' @param df A dataframe containing the measurement and grouping columns
#' @param sim A SIMPER result object as returned by `vegan::simper()`
#' @param x_axis Name of the x-axis grouping column (default: `'Year'`)
#' @param taxon_col Name of the taxon column (default: `'Taxon'`)
#' @param y_axis Name of the y-axis measurement column
#'   (default: `'Biovolume_per_mL'`)
#' @param agg_fun Aggregation function to apply per taxon per x group
#'   (default: `sum`)
#' @param p Significance threshold for filtering SIMPER taxa (default: `0.05`)
#' @param taxa_num Maximum number of top taxa to display (default: `6`)
#'
#' @return A `ggplot` object
#'
#' @importFrom dplyr filter arrange pull group_by summarise
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @importFrom tibble rownames_to_column
#' @importFrom stringr str_to_title
#' @export
plot_simper_summary <- function(df, sim,
                                x_axis = 'Year',
                                taxon_col = 'Taxon',
                                y_axis = 'Biovolume_per_mL',
                                agg_fun = sum,
                                p = 0.05,
                                taxa_num = 6) {
  taxon_sym <- sym(taxon_col)
  x_sym <- sym(x_axis)
  y_sym <- sym(y_axis)
  agg_name <- as.character(substitute(agg_fun))[[1]]

  top_taxa <- sim[[1]] %>%
    as.data.frame() %>%
    rownames_to_column(var = taxon_col) %>%
    arrange(desc(average)) %>%
    filter(p <= !!p) %>%
    slice_head(n = taxa_num) %>%
    pull(!!taxon_sym)

  df %>%
    filter((!!taxon_sym) %in% top_taxa) %>%
    group_by(!!taxon_sym, !!x_sym) %>%
    summarise(
      Value = agg_fun(!!y_sym, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    ggplot(aes(x = !!x_sym, y = Value, color = !!taxon_sym, group = !!taxon_sym)) +
    geom_line() +
    geom_point(size = 2) +
    labs(
      title = paste('Top SIMPER', y_axis),
      x = x_axis,
      y = paste(stringr::str_to_title(agg_name), 'of', as_label(y_sym)),
      color = 'Taxon'
    ) +
    theme_minimal()
}

#' @title Assign Delta regions to station coordinates
#'
#' @description
#' Converts a dataframe with `Latitude` and `Longitude` columns to a spatial
#' object, performs a spatial overlay against the Delta subregion geometry
#' (`delta_geom`), and returns an `sf` object with `Region` and `SubRegion`
#' columns added. Rows whose coordinates fall outside all regions are dropped.
#'
#' @param df A dataframe containing `Latitude` and `Longitude` columns
#' @param sf_delta Optional `sf` object of Delta subregion boundaries. If
#'   `NULL`, the package's built-in `delta_geom` dataset is used.
#'
#' @return An `sf` object with `Region` and `SubRegion` columns added and rows
#'   outside the Delta removed
#'
#' @importFrom sp CRS SpatialPointsDataFrame spTransform over
#' @importFrom dplyr filter
#' @importFrom tibble as_tibble
#' @export
convert_to_sf <- function(df, sf_delta = NULL) {

  # load in delta regional data
  if (is.null(sf_delta)) {
    if (!exists('delta_geom', envir = environment(), inherits = FALSE)) {
      utils::data('delta_geom', package = 'pespr', envir = environment())
    }
    sf_delta <- get('delta_geom', envir = environment())
  }

  # convert df to SpatialPointsDataFrame
  coords <- df[, c('Longitude', 'Latitude')]
  data   <- subset(df, select = -c(Latitude, Longitude))
  crs    <- CRS(SRS_string = 'EPSG:4326')

  spdf_wq <- SpatialPointsDataFrame(
    coords = coords,
    data   = data,
    proj4string = crs
  )

  # convert delta sf to Spatial
  spdf_delta <- as(sf_delta, 'Spatial')
  spdf_delta <- spTransform(spdf_delta, CRS('+init=epsg:4326 +proj=longlat'))

  # assign Region/SubRegion
  col_sr <- sp::over(spdf_wq, spdf_delta[, c('Region', 'SubRegion')])
  spdf_wq$Region <- col_sr$Region
  spdf_wq$SubRegion <- col_sr$SubRegion

  # convert to sf and tibble
  sf_wq <- st_as_sf(spdf_wq)
  sf_wq <- st_transform(sf_wq, st_crs(sf_delta))
  sf_wq <- sf_wq %>% filter(!is.na(Region))

  df_wq <- as_tibble(sf_wq)

  # return
  return(sf_wq)
}

#' @title Export community data in PRIMER format
#'
#' @description
#' Formats a community matrix and optional metadata into the layout expected by
#' PRIMER-e, then writes it to an Excel workbook. The output has the title in
#' the first row (merged across all columns), sample IDs (`R1`, `R2`, ...) in
#' the second row, community variables as rows below, and metadata factors
#' appended at the bottom after a blank row.
#'
#' @param df_meta A dataframe of sample metadata; one row per sample
#' @param df_com A community matrix dataframe or matrix; one row per sample,
#'   one column per taxon/variable. Must have the same number of rows as
#'   `df_meta`.
#' @param title Character string for the merged title cell
#'   (default: `'Environmental variables'`)
#' @param meta_as_factors Logical; if `TRUE`, metadata columns are appended
#'   below the community data as factor rows (default: `TRUE`)
#'
#' @return An `openxlsx` workbook object. Use `openxlsx::saveWorkbook()` to
#'   write it to disk.
#'
#' @importFrom dplyr mutate across
#' @importFrom tibble as_tibble rownames_to_column
#' @export
export_primer <- function(df_meta, df_com,
                          title = 'Environmental variables',
                          meta_as_factors = TRUE) {

  # check
  if (nrow(df_meta) != nrow(df_com))
    stop('df_meta and df_com must have same number of rows (samples).')

  # sample IDs
  sample_ids <- paste0('R', seq_len(nrow(df_meta)))

  # community matrix (variables as rows)
  df_com <- as_tibble(df_com)
  df_mat <- as.data.frame(t(df_com))
  colnames(df_mat) <- sample_ids
  df_mat <- rownames_to_column(df_mat, 'Variable')

  # row for sample IDs (R1, R2, …)
  row_sample_ids <- c('Variable', sample_ids)

  # optionally include metadata below
  df_factors <- NULL
  if (meta_as_factors) {
    df_factors <- df_meta %>%
      mutate(across(everything(), as.character)) %>%
      t() %>%
      as.data.frame(stringsAsFactors = FALSE)
    colnames(df_factors) <- sample_ids
    df_factors <- rownames_to_column(df_factors, 'Variable')
  }

  # prepare title + blank rows
  n_cols <- ncol(df_mat)
  row_title <- rep('', n_cols)
  row_title[1] <- title
  row_blank <- rep('', n_cols)

  # assemble in correct PRIMER order:
  # title, sample IDs, community data, blank, metadata
  df_primer <- rbind(
    row_title,
    row_sample_ids,
    df_mat,
    row_blank,
    df_factors
  )

  # replace header cell with NA as PRIMER expects
  df_primer[2, 1] <- NA

  # create workbook
  wb <- createWorkbook()
  addWorksheet(wb, 'Sheet1')

  # write to Excel
  writeData(wb, 'Sheet1', df_primer, colNames = FALSE)

  # merge title row
  mergeCells(wb, 'Sheet1', cols = 1:n_cols, rows = 1)

  return(wb)
}


# Plot Cells >= Units -----------------------------------------------------

#' @title Plot raw cell count vs unit abundance differences over time
#'
#' @description
#' For a specified taxon, computes the difference between total cell count and
#' unit abundance columns and plots it as a scatter plot over time. Points where
#' cells are less than units are highlighted in red. Shaded regions and dashed
#' vertical lines mark known periods of data quality issues.
#'
#' @param df_data A dataframe containing `Date`, and the specified cell and
#'   unit columns
#' @param taxon Taxon name to use in the plot title (default:
#'   `'Microcystis aeruginosa'`)
#' @param cell_col Name of the total cell count column
#'   (default: `'Total Number of Cells'`)
#' @param unit_col Name of the unit abundance column
#'   (default: `'Unit Abundance (# of Natural Units)'`)
#' @param x_range Optional two-element character vector of dates (`'YYYY-MM-DD'`)
#'   to zoom the x-axis via `coord_cartesian()`
#' @param y_range Optional two-element numeric vector to zoom the y-axis
#'
#' @return A `ggplot` object
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_rect geom_point scale_color_manual geom_hline geom_vline labs theme_minimal coord_cartesian
#' @export
plot_cellunit_check <- function(df_data,
                           taxon = 'Microcystis aeruginosa',
                           cell_col = 'Total Number of Cells',
                           unit_col = 'Unit Abundance (# of Natural Units)',
                           x_range = NULL,
                           y_range = NULL) {

  # convert string column names to symbols
  cell_sym <- sym(cell_col)
  unit_sym <- sym(unit_col)

  # filter and compute difference
  df_micro <- df_data %>%
    mutate(Diff = !!cell_sym - !!unit_sym)

  # define shaded regions
  df_shade <- data.frame(
    xmin = as.Date(c('2013-09-01', '2013-11-01', '2014-01-01', '2021-03-01')),
    xmax = as.Date(c('2013-09-30', '2013-11-30', '2020-12-31', '2021-10-31'))
  )

  # define vertical lines
  vlines <- as.Date(c('2013-09-01','2013-09-30','2013-11-01','2013-11-30','2014-01-01','2020-12-31','2021-03-01','2021-10-31'))

  # build plot
  plt <- ggplot(df_micro, aes(x = Date, y = Diff)) +
    geom_rect(data = df_shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = 'gray85', alpha = 0.5) +
    geom_point(aes(color = Diff < 0)) +
    scale_color_manual(values = c('FALSE' = 'steelblue', 'TRUE' = 'red'), guide = 'none') +
    geom_hline(yintercept = 0, linetype = 'solid', color = 'black') +
    geom_vline(xintercept = vlines, linetype = 'dashed', color = 'gray40') +
    labs(title = paste('Differences for',taxon), x = 'Date', y = 'TotalCells - UnitAbundance') +
    theme_minimal()

  # apply zoom
  if (!is.null(x_range) || !is.null(y_range)) {
    plt <- plt + coord_cartesian(xlim = as.Date(x_range), ylim = y_range)
  }

  return(plt)
}
