library(miRNARefine)

test_that("adaptiveFiltering input type check", {
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
