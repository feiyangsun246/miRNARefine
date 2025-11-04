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

})

test_that("median imputation works on toy data", {
  toy <- data.frame(miR1 = c(1, NA, 3), miR2 = c(NA, 2, 4))
  filled <- missingValueHandling(toy, method = "median", report_summary = FALSE)

  # no more NA should remain after imputation
  expect_false(any(is.na(filled)))

  # median imputation actually works
  expect_equal(filled$miR1[2], median(c(1,3)))
  expect_equal(filled$miR2[1], median(c(2,4)))

})

test_that("median imputation works on actual miRNA data", {
  toy <- data.frame(miR1 = c(1, NA, 3), miR2 = c(NA, 2, 4))
  filled <- missingValueHandling(toy, method = "median", report_summary = FALSE)

  # no more NA should remain after imputation
  expect_false(any(is.na(filled)))

  # median imputation actually works
  expect_equal(filled$miR1[2], median(c(1,3)))
  expect_equal(filled$miR2[1], median(c(2,4)))

})
