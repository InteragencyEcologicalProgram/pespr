
# Plot NMDS ---------------------------------------------------------------

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
  
  if (compute_dist) {
    dist_obj <- vegdist(m_com, method = distance)
    print('here')
  }
  
  
  return_list <- list(
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
  
  return(return_list)
}

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

# map regions
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

# export primer
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