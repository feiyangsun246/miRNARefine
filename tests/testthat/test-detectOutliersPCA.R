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
  expect_true(is.list(result))
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


test_that("Check if detectOutliersPCA correctly identifies outliers on
          simple data", {

  # Clear sample outliers
  set.seed(624)
  df <- data.frame(
    miR1 = c(rnorm(9, 0, 1), 100),
    miR2 = c(rnorm(9, 0, 1), 100),
    miR3 = c(rnorm(9, 0, 1), 100),
    miR4 = c(rnorm(9, 0, 1), 100)
  )
  result <- detectOutliersPCA(df, z_threshold = 2, report_summary = FALSE)
  expect_true(is.list(result))
  s_outliers <- which(result$sample_outliers)

  # Whether 10th row is detected as an outlier
  expect_true(10 %in% s_outliers)

  # Whether at least one miRNA in the 10th row is flagged
  expect_true(any(result$miRNA_outliers[10, ]))

  # Whether other rows are flagged
  expect_false(any(result$sample_outliers[1:9]))
  expect_false(any(result$miRNA_outliers[1:9, ]))

})


test_that("Check if detectOutliersPCA stops when NA values present", {
  df <- data.frame(
    x1 = c(1, 2, NA, 4),
    x2 = c(1, 2, 3, 4)
  )
  expect_error(detectOutliersPCA(df, report_summary = FALSE),
               "Dataset contains missing values")
})

