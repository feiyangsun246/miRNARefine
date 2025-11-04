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

  # Empty input
  df_empty <- data.frame()
  expect_error(adaptiveFiltering(df_empty),
               "Empty dataframe input")
})


test_that("Check if detectOutliersPCA output structure is correct", {

  result <- detectOutliersPCA(miRNASeq1,
                              report_summary = FALSE)

  # Whether result contain such elements
  expect_true("pca_res" %in% names(result))
  expect_true("miRNA_outliers" %in% names(result))
  expect_true(is.logical(result$miRNA_outliers))
})


test_that("Check if detectOutliersPCA handles singular covariance matrix
          correctly", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    y = c(2, 4, 6, 8)  # singular
  )

  expect_warning(result <- detectOutliersPCA(df, report_summary = FALSE),
                 "Covariance matrix is singular")

  expect_true(all(c("Sample", "IsOutlier") %in% colnames(result)))

  # Should all return FALSE
  expect_false(any(result$IsOutlier))
})


