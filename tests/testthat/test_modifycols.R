
# remove_nodata -----------------------------------------------------------

test_that('remove_nodata works correctly', {
  # 1. tibble with valid data (no zeros or all-NA rows)
  df1 <- tibble(
    Cells_per_mL = c(1, 2, 3),
    Units_per_mL = c(4, 5, 6),
    Biovolume_per_mL = c(7, 8, 9)
  )
  res1 <- remove_nodata(df1)
  expect_equal(nrow(res1), 3)
  expect_null(attr(res1, 'log')$nodata)
  
  # 2. tibble where all are NA (should remove all rows)
  df2 <- tibble(
    Cells_per_mL = c(NA, NA),
    Units_per_mL = c(NA, NA),
    Biovolume_per_mL = c(NA, NA)
  )
  res2 <- remove_nodata(df2)
  expect_equal(nrow(res2), 0)
  expect_equal(nrow(attr(res2, 'log')$nodata), 2)
  
  # 3. tibble where some are NA (keep partial rows)
  df3 <- tibble(
    Cells_per_mL = c(NA, 2, NA),
    Units_per_mL = c(3, NA, 4),
    Biovolume_per_mL = c(NA, 5, 6)
  )
  res3 <- remove_nodata(df3)
  expect_equal(nrow(res3), 3)  # none removed since no row has all NA
  expect_equal(sum(is.na(res3$Cells_per_mL)), 2)  # still keeps NAs
  expect_null(attr(res3, 'log')$nodata)
  
  # 4. tibble where some are 0 (should become NA but not removed unless all NA)
  df4 <- tibble(
    Cells_per_mL = c(0, 1),
    Units_per_mL = c(2, 0),
    Biovolume_per_mL = c(3, 4)
  )
  res4 <- remove_nodata(df4)
  expect_true(all(is.na(res4$Cells_per_mL[1])))
  expect_true(all(is.na(res4$Units_per_mL[2])))
  expect_equal(nrow(res4), 2)
  
  # 5. tibble where some are 0.0 (should behave same as 0)
  df5 <- tibble(
    Cells_per_mL = c(0.0, 1.1),
    Units_per_mL = c(2.2, 0.0),
    Biovolume_per_mL = c(3.3, 4.4)
  )
  res5 <- remove_nodata(df5)
  expect_true(is.na(res5$Cells_per_mL[1]))
  expect_true(is.na(res5$Units_per_mL[2]))
  expect_equal(nrow(res5), 2)
  
  # 6. tibble where some are 0.00 (should behave same as 0)
  df6 <- tibble(
    Cells_per_mL = c(0.00, 1.1),
    Units_per_mL = c(2.2, 0.00),
    Biovolume_per_mL = c(3.3, 0.00)
  )
  res6 <- remove_nodata(df6)
  expect_true(is.na(res6$Cells_per_mL[1]))
  expect_true(is.na(res6$Units_per_mL[2]))
  expect_true(is.na(res6$Biovolume_per_mL[2]))
  expect_equal(nrow(res6), 2)
})
