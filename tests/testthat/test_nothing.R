library(testthat)
library(FeralCatEradication)

test_that("all is ok", {
  expected_value = 1
  obtained_value = return_one()
  expect_equal(expected_value, obtained_value)
})
