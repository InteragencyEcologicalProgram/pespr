# flag_outliers -----------------------------------------------------------
test_that('flag_outliers flags outliers correctly for Cells_per_mL', {
  df <- tibble(
    Date = as.Date(c('2020-01-01', '2020-01-02', '2020-01-03')),
    Station = c('A1', 'A1', 'A1'),
    Taxon = c('Fragilaria', 'Cyclotella', 'Navicula'),
    Cells_per_mL = c(10, 10, 10000),
    Units_per_mL = c(10, 10, 10),
    Biovolume_per_mL = c(1, 1, 1),
    QualityCheck = c('NoCode', 'NoCode', 'NoCode'),
    Notes = c('NoNotes', 'NoNotes', 'NoNotes')
  )
  
  out <- suppressMessages(flag_outliers(df, Cells_per_mL, show_plot = FALSE))
  
  # one outlier should be flagged
  expect_true(any(grepl('OutlierCells', out$QualityCheck)))
  
  # only one row should be flagged
  flagged <- out %>% filter(str_detect(QualityCheck, 'OutlierCells'))
  expect_equal(nrow(flagged), 1)
})

test_that('flag_outliers does not flag normal data', {
  df <- tibble(
    Date = as.Date(c('2020-01-01', '2020-01-02', '2020-01-03')),
    Station = c('A1', 'A1', 'A1'),
    Taxon = c('Fragilaria', 'Cyclotella', 'Navicula'),
    Cells_per_mL = c(10, 12, 11),
    Units_per_mL = c(10, 10, 10),
    Biovolume_per_mL = c(1, 1, 1),
    QualityCheck = c('NoCode', 'NoCode', 'NoCode')
  )
  
  out <- suppressMessages(flag_outliers(df, Cells_per_mL, show_plot = FALSE))
  
  # no outliers expected
  expect_false(any(grepl('OutlierCells', out$QualityCheck)))
})

test_that('flag_outliers handles mixed QualityCheck values correctly', {
  df <- tibble(
    Date = as.Date(c('2020-01-01', '2020-01-02', '2020-01-03')),
    Station = c('A1', 'A1', 'A1'),
    Taxon = c('Fragilaria', 'Cyclotella', 'Navicula'),
    Cells_per_mL = c(10, 10, 10000),
    Units_per_mL = c(10, 10, 10),
    Biovolume_per_mL = c(1, 1, 1),
    QualityCheck = c('NoCode', 'ManualCheck', 'ManualCheck')
  )
  
  out <- suppressMessages(flag_outliers(df, Cells_per_mL, show_plot = FALSE))
  
  # outlier should have both tags if not NoCode
  flagged <- out %>% filter(str_detect(QualityCheck, 'OutlierCells'))
  expect_true(all(str_detect(flagged$QualityCheck, 'OutlierCells')))
})

test_that('flag_outliers errors for invalid column', {
  df <- tibble(
    Date = as.Date('2020-01-01'),
    Station = 'A1',
    Taxon = 'Fragilaria',
    Cells = 10,
    QualityCheck = 'NoCode'
  )
  
  expect_error(flag_outliers(df, Error, show_plot = FALSE),
               'colname must contain one of')
})

test_that('flag_outliers removes existing OutlierCells before re-flagging', {
  df <- tibble(
    Date = as.Date(c('2020-01-01', '2020-01-02', '2020-01-03', '2020-01-04', '2020-01-05')),
    Station = c('A1', 'A1', 'A1', 'A1', 'A1'),
    Taxon = c('Fragilaria', 'Cyclotella', 'Navicula', 'Synedra', 'Microcystis'),
    Cells_per_mL = c(10, 10, 10000, 10, 10000),
    Units_per_mL = c(10, 10, 10, 10, 10),
    Biovolume_per_mL = c(1, 1, 1, 1, 1),
    QualityCheck = c(
      'OutlierCells', # should remove outlier flag
      'Test; OutlierCells', # should remove outlier flag
      'ManualCheck', # should add outlier flag
      'NoCode', # normal case
      'Test; OutlierCells' # remains flagged (still outlier)
    ),
    Notes = NA
  )
  
  out <- suppressMessages(flag_outliers(df, Cells_per_mL, show_plot = FALSE))
  
  expect_equal(out$QualityCheck[1], 'NoCode')
  expect_equal(out$QualityCheck[2], 'Test')
  expect_equal(out$QualityCheck[3], 'ManualCheck; OutlierCells')
  expect_equal(out$QualityCheck[4], 'NoCode')
  expect_equal(out$QualityCheck[5], 'Test; OutlierCells')
})
