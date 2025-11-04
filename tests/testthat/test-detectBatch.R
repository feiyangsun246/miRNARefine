library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if detectBatch error upon invalid dataset input", {

  # Handling missing values to avoid unexpected errors
  filled1 <- missingValueHandling(miRNAdata = miRNASeq1, method = "median",
                                  report_summary = FALSE)
  filled2 <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
                                  report_summary = FALSE)

  # Valid input
  expect_silent(compareNormalization(miRNAdata = filled2,
                                     report_summary = FALSE))

  # Invalid input: Vector
  expect_error(detectBatch(miRNAdata = miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(detectBatch(miRNAdata = as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(detectBatch(miRNAdata = 624),
               "`miRNAdata` must be a data frame or matrix")

  # Empty input
  df_empty <- data.frame()
  expect_error(compareNormalization(miRNAdata = df_empty),
               "Empty dataframe input")
})
