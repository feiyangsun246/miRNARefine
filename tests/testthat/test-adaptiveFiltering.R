library(miRNARefine)

test_that("Check if input is valid", {
  data("miRNASeq1")
  data("miRNASeq2")

  # Valid input
  expect_silent(adaptiveFiltering(miRNAdata = miRNASeq1,
                                  report_summary = FALSE))

  # Valid input
  expect_silent(adaptiveFiltering(miRNAdata = miRNASeq2,
                                  report_summary = FALSE))

  # Invalid input: Vector
  expect_error(adaptiveFiltering(miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(adaptiveFiltering(as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(adaptiveFiltering(42),
               "`miRNAdata` must be a data frame or matrix")

})

test_that("Check if output is a data.frame with expected dimensions", {
  data("miRNASeq1")
  data("miRNASeq2")
  filtered1 <- adaptiveFiltering(miRNASeq1)

  # Output class check for miRNASeq1
  expect_s3_class(filtered1, "data.frame")

  # Number of columns and rows in output for miRNASeq1
  expect_true(nrow(filtered1) == nrow(miRNASeq1))
  expect_true(ncol(filtered1) <= ncol(miRNASeq1))

  filtered2 <- adaptiveFiltering(miRNASeq2)

  # Output class check for miRNASeq2
  expect_s3_class(filtered2, "data.frame")

  # Number of columns and rows in output for miRNASeq2
  expect_true(nrow(filtered2) == nrow(miRNASeq2))
  expect_true(ncol(filtered2) <= ncol(miRNASeq2))


})
