library(miRNARefine)

test_that("Check if input is valid", {
  data("miRNASeq1")
  data("miRNASeq2")

  # Valid input
  expect_silent(detectOutliersPCA(miRNAdata = miRNASeq1,
                                  report_summary = FALSE))

  # Valid input
  expect_silent(detectOutliersPCA(miRNAdata = miRNASeq2,
                                  report_summary = FALSE))

  # Invalid input: Vector
  expect_error(detectOutliersPCA(miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(detectOutliersPCA(as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(detectOutliersPCA(42),
               "`miRNAdata` must be a data frame or matrix")

})
