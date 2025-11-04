library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if missingValueHandling error upon invalid user input", {

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

  # Empty input
  df_empty <- data.frame()
  expect_error(adaptiveFiltering(df_empty),
               "Empty dataframe input")
})


test_that("factor columns are converted to numeric", {

  toy <- data.frame(miR1 = factor(c(1, 2, NA)))
  filled <- missingValueHandling(toy, method="median",
                                 report_summary = FALSE)
  # Whether all NA handled
  expect_false(any(is.na(filled)))

  # Whether conversion is completed
  expect_true(is.numeric(filled$miR1))
})


test_that("Check if median imputation works on toy data", {

  toy <- data.frame(miR1 = c(1, NA, 3), miR2 = c(NA, 2, 4))
  filled <- missingValueHandling(toy, method = "median", report_summary = FALSE)

  # No more NA should remain after imputation
  expect_false(any(is.na(filled)))

  # Whether median imputation actually works
  expect_equal(filled$miR1[2], median(c(1,3)))
  expect_equal(filled$miR2[1], median(c(2,4)))
})


test_that("Check if median imputation works on actual miRNA data", {

  # Test on a subpart of miRNASeq1
  subpart <- miRNASeq1[1:10, 1:2]
  subpart[1, 1] <- NA

  # Whether the value calculated is as expected
  filled_sub <- missingValueHandling(subpart, method = "median",
                                     report_summary = FALSE)
  expect_equal(filled_sub[1, 1],
               median(filled_sub[2:10, 1]))

  # Test on the whole set of miRNASeq1
  filled <- missingValueHandling(miRNASeq1, method = "median",
                                 report_summary = FALSE)

  # Whether all missing data filled
  expect_false(any(is.na(filled)))

  # Whether dataset structure remains the same
  expect_equal(ncol(filled), ncol(miRNASeq1))
  expect_equal(nrow(filled), nrow(miRNASeq1))
})


test_that("Check if mean imputation works on toy data", {

  toy <- data.frame(miR1 = c(1, NA, 3))
  filled <- missingValueHandling(toy, method = "mean",
                                 report_summary = FALSE)

  # Whether mean imputation actually works
  expect_equal(filled$miR1[2], mean(c(1,3)))
})


test_that("Check if mean imputation works on actual miRNA data", {

  # Test on a subpart of miRNASeq1
  subpart <- miRNASeq1[1:10, 1:2]
  subpart[1, 1] <- NA

  # Whether the value calculated is as expected
  filled_sub <- missingValueHandling(subpart, method = "mean",
                                     report_summary = FALSE)
  expect_equal(filled_sub[1, 1],
               mean(filled_sub[2:10, 1]))

  # Test on the whole set of miRNASeq1
  filled <- missingValueHandling(miRNASeq1, method = "mean",
                                 report_summary = FALSE)

  # Whether all missing data filled
  expect_false(any(is.na(filled)))

  # Whether dataset structure remains the same
  expect_equal(ncol(filled), ncol(miRNASeq1))
  expect_equal(nrow(filled), nrow(miRNASeq1))
})


test_that("Check if knn imputation works", {

  # Whether package impute is installed
  skip_if_not_installed("impute")

  toy <- matrix(c(1, NA, 3, 4, 3, 2), nrow=3, ncol=2)
  filled <- missingValueHandling(toy, method="knn", k=2,
                                 report_summary = FALSE)

  # Whether no NA remains
  expect_false(any(is.na(filled)))

  # Whether structure of dataset remains the same
  expect_equal(dim(filled), dim(toy))
})

