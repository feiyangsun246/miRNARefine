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
  batch1 <- rep(1:3, each = 10)
  expect_silent(detectBatch(miRNAdata = filled2,
                            batch = batch1,
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
  expect_error(detectBatch(miRNAdata = df_empty),
               "Empty dataframe input")
})


test_that("Check if error appears when length of batch doesn't match number
          of samples", {

  # Invalid batch input: vector
  expect_error(detectBatch(miRNAdata = matrix(1:10, 5, 2), batch = 1:3,
                           report_summary = FALSE),
               "Length of batch")

  # Invalid batch input: a single number
  expect_error(detectBatch(miRNAdata = matrix(1:10, 5, 2), batch = 1,
                           report_summary = FALSE),
               "Length of batch")

  # Valid batch input: NULL
  expect_silent(detectBatch(miRNAdata = matrix(1:10, 5, 2), batch = NULL,
                           report_summary = FALSE))
})


test_that("Check if situation is correctly handled when batch factor is
          provided as NULL", {

  filled2 <- missingValueHandling(miRNAdata = miRNASeq2,
                                  method = "median",
                                  report_summary = FALSE)

  expect_warning(res <- detectBatch(miRNAdata = filled2,
                                    correct = TRUE,
                                    report_summary = FALSE))
  expect_equal(res$corrected_data, filled2)
})


test_that("Check if output structure is correct", {

  set.seed(624)
  df <- matrix(rnorm(20), nrow = 4, ncol = 5)
  batch1 <- rep(1:2, each = 2)
  res <- detectBatch(miRNAdata = df, batch = batch1,
                     correct = TRUE, report_summary = FALSE)

  # Whether result is a list
  expect_type(res, "list")

  # Whether result contains certain domains
  expect_true(all(c("corrected_data", "batch_effects") %in% names(res)))
  expect_equal(dim(res$corrected), dim(df))
})


test_that("Check if Batch effect correction works on simple datasets", {

  set.seed(624)
  df <- matrix(rnorm(20), 5, 4)
  bat <- c(1,1,2,2,2)

  # Add batch effect
  df[bat == 2, ] <- df[bat == 2, ] + 5

  res <- detectBatch(miRNAdata = df, batch = bat, correct = TRUE,
                     report_summary = FALSE)

  # Before correction, mean difference between batches
  mean_diff_before <- mean(colMeans(df[bat == 1, , drop=FALSE])) -
                      mean(colMeans(df[bat == 2, , drop=FALSE]))

  # After correction
  batch1_rows <- res$corrected_data[bat == 1, , drop=FALSE]
  batch2_rows <- res$corrected_data[bat == 2, , drop=FALSE]

  mean_diff_after <- mean(colMeans(batch1_rows)) - mean(colMeans(batch2_rows))

  expect_lt(abs(mean_diff_after), abs(mean_diff_before))

  expect_type(res, "list")
  expect_true(all(c("corrected_data", "batch_effects") %in% names(res)))
  expect_equal(dim(res$corrected), dim(df))
})


test_that("Check if Batch effect correction works on real miRNA datasets", {
  bat <- rep(1:2, each = 5)
  filled1 <- missingValueHandling(miRNAdata = miRNASeq1, method = "median",
                                  report_summary = FALSE)

  filled1 <- filled1[, colSums(filled1 != 0) > 0]

  # Add batch effect
  filled1[bat == 2, ] <- filled1[bat == 2, ] + 5

  res <- suppressMessages(detectBatch(miRNAdata = filled1,
                                      batch = bat,
                                      correct = TRUE,
                                      report_summary = FALSE))

  # Before correction, mean difference between batches
  mean_diff_before <- mean(colMeans(filled1[bat == 1, , drop=FALSE])) -
                      mean(colMeans(filled1[bat == 2, , drop=FALSE]))

  # After correction
  batch1_rows <- res$corrected_data[bat == 1, , drop=FALSE]
  batch2_rows <- res$corrected_data[bat == 2, , drop=FALSE]

  mean_diff_after <- mean(colMeans(batch1_rows)) - mean(colMeans(batch2_rows))

  expect_lt(abs(mean_diff_after), abs(mean_diff_before))

  expect_type(res, "list")
  expect_true(all(c("corrected_data", "batch_effects") %in% names(res)))
  expect_equal(dim(res$corrected), dim(filled1))
})

