library(miRNARefine)

test_that("Check if missingValueHandling error upon invalid user input", {
  data("miRNASeq1")
  data("miRNASeq2")

  # Valid input
  expect_silent(missingValueHandling(miRNAdata = miRNASeq1,
                                  report_summary = FALSE))

  # Valid input
  expect_silent(missingValueHandling(miRNAdata = miRNASeq2,
                                  report_summary = FALSE))

  # Invalid input: Vector
  expect_error(missingValueHandling(miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(missingValueHandling(as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(missingValueHandling(624),
               "`miRNAdata` must be a data frame or matrix")

})
