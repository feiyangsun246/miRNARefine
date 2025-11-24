library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if adaptiveFiltering error upon invalid dataset input", {

  # Valid input
  expect_silent(adaptiveFiltering(miRNAdata = miRNASeq1,
                                  report_summary = FALSE))

  # Valid input
  expect_silent(adaptiveFiltering(miRNAdata = miRNASeq2,
                                  report_summary = FALSE))

  # Invalid input: Vector
  expect_error(adaptiveFiltering(miRNAdata = miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(adaptiveFiltering(miRNAdata = as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(adaptiveFiltering(miRNAdata = 624),
               "`miRNAdata` must be a data frame or matrix")

  # Empty input
  df_empty <- data.frame()
  expect_error(adaptiveFiltering(miRNAdata = df_empty),
               "Empty dataframe input")
})

test_that("Check if adaptiveFiltering error upon negative thresholds", {

  # Use dataset provided in the package
  expect_silent(adaptiveFiltering(miRNAdata = miRNASeq1,
                                  report_summary = FALSE))
  # Negative min_expression
  expect_error(adaptiveFiltering(miRNAdata = miRNASeq1,
                                 min_expression = -200,
                                 report_summary = FALSE),
               "thresholds must all be non-negative.")

  # Negative min_variance
  expect_error(adaptiveFiltering(miRNAdata = miRNASeq1,
                                 min_variance = -200,
                                 report_summary = FALSE),
               "thresholds must all be non-negative.")

  # Negative max_na
  expect_error(adaptiveFiltering(miRNAdata = miRNASeq1,
                                 max_na = -200,
                                 report_summary = FALSE),
               "thresholds must all be non-negative.")
})


test_that("Check if output of adaptiveFiltering is a data.frame with expected
          dimensions", {

  data("miRNASeq1")
  data("miRNASeq2")
  filtered1 <- adaptiveFiltering(miRNAdata = miRNASeq1, report_summary = FALSE)

  # Output class check for miRNASeq1
  expect_s3_class(filtered1, "data.frame")

  # Number of columns and rows in output for miRNASeq1
  expect_true(nrow(filtered1) == nrow(miRNASeq1))
  expect_true(ncol(filtered1) <= ncol(miRNASeq1))

  filtered2 <- adaptiveFiltering(miRNAdata = miRNASeq2, report_summary = FALSE)

  # Output class check for miRNASeq2
  expect_s3_class(filtered2, "data.frame")

  # Number of columns and rows in output for miRNASeq2
  expect_true(nrow(filtered2) == nrow(miRNASeq2))
  expect_true(ncol(filtered2) <= ncol(miRNASeq2))
})


test_that("Check if low-expression and low-variance miRNAs are filtered by
          adaptiveFiltering", {

  # Processing with toy data
  toy <- data.frame(
    miR1 = c(0.1, 0.2, 0.1, 0.2),
    miR2 = c(1, 1, 1, 1),
    miR3 = c(1, 2, 3, 4)
  )
  result <- adaptiveFiltering(miRNAdata = toy,
                              min_expression = 0.3,
                              min_variance = 0.1,
                              max_na = 0,
                              report_summary = FALSE)

  # Check whether miR1 and miR2 are removed
  expect_true("miR3" %in% colnames(result))
  expect_false("miR1" %in% colnames(result))
  expect_false("miR2" %in% colnames(result))
})


test_that("Check if adaptiveFiltering removes columns with too many
          missing values", {

  # Processing with toy data
  toy <- data.frame(
    miR1 = c(1, 2, 3, NA),      # 25% NA
    miR2 = c(NA, NA, NA, 1)     # 75% NA
  )

  # Check whether miR2 is removed
  result <- adaptiveFiltering(miRNAdata = toy,
                              max_na = 0.5,
                              report_summary = FALSE)
  expect_true("miR1" %in% colnames(result))
  expect_false("miR2" %in% colnames(result))
})

