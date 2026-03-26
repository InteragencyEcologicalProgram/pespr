utils::globalVariables(c(
  # operators / tidy eval
  '%>%', ':=', '.', '.data',
  # data.table internals
  '.N', '.SD',
  # dplyr NSE column names
  'AlgalGroup', 'Cells_per_mL', 'Class', 'CurrentTaxon', 'Date',
  'Date_new', 'Date_orig', 'Debris', 'DepthType', 'Diff', 'Dimension',
  'Ending Date', 'Event_ID', 'GALD', 'GALD 1', 'Genus', 'Kingdom',
  'Latitude', 'Longitude', 'Month', 'MonthYear', 'N', 'NMDS1', 'NMDS2',
  'Notes', 'NumStations', 'OldTaxon', 'OrigTaxon', 'Phylum', 'PhytoForm',
  'PureTaxon', 'PureTaxon_check', 'QC_3', 'QC_4', 'QC_5', 'QualityCheck',
  'QualityCheck_new', 'Region', 'SampleDate', 'SampleDepth', 'SampleTime',
  'Sampling Method', 'Species', 'Starting Date', 'Station', 'StationCode',
  'Survey', 'Taxon', 'Time', 'Time_chr', 'Time_fixed', 'Time_new',
  'Time_orig', 'TimesPerYear', 'Units_per_mL', 'Value', 'Year',
  'Colony/Filament/Individual Group Code',
  # internal computed variables
  '.orig_taxon', '.paren_content', '.time_hms', 'average', 'both_4_5',
  'datetime_local', 'datetime_pst', 'dep', 'depth_code', 'depth_num',
  'end', 'grab_code', 'grab_num', 'has_class', 'has_king', 'has_phyl',
  'log_val', 'n', 'n_distinct_values', 'n_effort', 'outlier',
  'survey_short', 'time_missing', 'total', 'value', 'x_index', 'xmax',
  'xmin', 'z_robust', 'start',
  # package data objects
  'phyto_taxonomy', 'phyto_typos', 'survey_metadata',
  # base R functions that generate false positives
  'as', 'setNames'
))
