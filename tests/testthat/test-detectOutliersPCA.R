library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if detectOutliersPCA error upon invalid user input", {

  # Handling missing values to avoid unexpected errors
  filled1 <- missingValueHandling(miRNASeq1, method = "median",
                                  report_summary = FALSE)
  filled2 <- missingValueHandling(miRNASeq2, method = "median",
                                  report_summary = FALSE)

  # Valid input
  expect_silent(detectOutliersPCA(miRNAdata = filled1,
                                  report_summary = FALSE))

  # Valid input
  expect_silent(detectOutliersPCA(miRNAdata = filled2,
                                  report_summary = FALSE))

  # Invalid input: Vector
  expect_error(detectOutliersPCA(filled1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(detectOutliersPCA(as.list(filled1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(detectOutliersPCA(624),
               "`miRNAdata` must be a data frame or matrix")

})

