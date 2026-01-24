
# remove_nodata -----------------------------------------------------------

test_that('remove_nodata keeps valid rows unchanged', {
  df <- tibble(
    Taxon = c('Fragilaria', 'Fragilaria', 'Fragilaria'),
    Cells_per_mL = c(1, 2, 3),
    Units_per_mL = c(4, 5, 6),
    Biovolume_per_mL = c(7, 8, 9)
  )
  res <- remove_nodata(df)
  expect_equal(nrow(res), 3)
  expect_null(attr(res, 'log')$nodata)
})

test_that('remove_nodata removes rows where all density cols are NA', {
  df <- tibble(
    Taxon = c('Fragilaria', 'Fragilaria'),
    Cells_per_mL = c(NA, NA),
    Units_per_mL = c(NA, NA),
    Biovolume_per_mL = c(NA, NA)
  )
  res <- remove_nodata(df)
  expect_equal(nrow(res), 0)
  expect_equal(nrow(attr(res, 'log')$nodata), 2)
})

test_that('remove_nodata keeps partially missing rows', {
  df <- tibble(
    Taxon = c('Fragilaria', 'Fragilaria', 'Fragilaria'),
    Cells_per_mL = c(NA, 2, NA),
    Units_per_mL = c(3, NA, 4),
    Biovolume_per_mL = c(NA, 5, 6)
  )
  res <- remove_nodata(df)
  expect_equal(nrow(res), 3)  # none removed since no row has all NA
  expect_equal(sum(is.na(res$Cells_per_mL)), 2)  # still keeps NAs
  expect_null(attr(res, 'log')$nodata)
})

test_that('remove_nodata converts zero values to NA', {
  df <- tibble(
    Taxon = c('Fragilaria', 'Fragilaria'),
    Cells_per_mL = c(0, 1),
    Units_per_mL = c(2, 0),
    Biovolume_per_mL = c(3, 4)
  )
  res <- remove_nodata(df)
  expect_true(all(is.na(res$Cells_per_mL[1])))
  expect_true(all(is.na(res$Units_per_mL[2])))
  expect_equal(nrow(res), 2)
})

test_that('remove_nodata treats 0.0 and 0.00 as zeros', {
  df <- tibble(
    Taxon = c('Fragilaria', 'Fragilaria'),
    Cells_per_mL = c(0.00, 1.1),
    Units_per_mL = c(2.2, 0.0),
    Biovolume_per_mL = c(3.3, 0.00)
  )
  res <- remove_nodata(df)
  expect_true(is.na(res$Cells_per_mL[1]))
  expect_true(is.na(res$Units_per_mL[2]))
  expect_true(is.na(res$Biovolume_per_mL[2]))
  expect_equal(nrow(res), 2)
})
