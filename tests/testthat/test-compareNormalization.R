library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if compareNormalization error upon invalid dataset input", {

  # Handling missing values to avoid unexpected errors
  filled1 <- missingValueHandling(miRNAdata = miRNASeq1, method = "median",
                                  report_summary = FALSE)
  filled2 <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
                                  report_summary = FALSE)

  # Valid input
  expect_silent(compareNormalization(miRNAdata = filled2,
                                     report_summary = FALSE))

  # Invalid input: Vector
  expect_error(compareNormalization(miRNAdata = miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(compareNormalization(miRNAdata = as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(compareNormalization(miRNAdata = 624),
               "`miRNAdata` must be a data frame or matrix")

  # Empty input
  df_empty <- data.frame()
  expect_error(compareNormalization(miRNAdata = df_empty),
               "Empty dataframe input")
})


test_that("Check if compareNormalization stops when NA values present", {

  df <- data.frame(
    x1 = c(1, 2, NA, 4),
    x2 = c(1, 2, 3, 4)
  )
  expect_error(compareNormalization(miRNAdata = df,
                                    report_summary = FALSE),
               "Dataset contains missing values")
})


test_that("Check if compareNormalization error upon invalid methods input", {

  filled2 <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
                                  report_summary = FALSE)
  # Invalid method: empty input
  expect_error(compareNormalization(miRNAdata = filled2,
                                    methods = "",
                                    report_summary = FALSE),
               "`methods` must be one or more of 'log2', 'zscore', or 'quantile'")

  # Invalid method: something other than the three options
  expect_error(compareNormalization(miRNAdata = filled2,
                                    methods = "something",
                                    report_summary = FALSE),
               "`methods` must be one or more of 'log2', 'zscore', or 'quantile'")

  # Invalid method: invalid combinations
  expect_error(compareNormalization(miRNAdata = filled2,
                                    methods = c("log2", "something"),
                                    report_summary = FALSE),
               "`methods` must be one or more of 'log2', 'zscore', or 'quantile'")

  # Invalid method: numeric
  expect_error(compareNormalization(miRNAdata = filled2,
                                    methods = 624,
                                    report_summary = FALSE),
               "`methods` must be one or more of 'log2', 'zscore', or 'quantile'")
})


test_that("Check if output structure is correct", {

  set.seed(624)
  df <- matrix(abs(rnorm(20)), nrow = 5, ncol = 4)
  res <- compareNormalization(miRNAdata = df,
                              report_summary = FALSE)

  # Whether output is a list
  expect_type(res, "list")

  # Whether output contains certain domains
  expect_true(all(c("normalized", "best_method") %in% names(res)))
  expect_true(all(c("log2", "zscore", "quantile") %in% names(res$normalized)))

  expect_equal(dim(res$normalized$log2), dim(df))
})


test_that("Check if best method selection", {

  set.seed(624)
  df <- matrix(abs(rnorm(20)), nrow = 5)
  res <- compareNormalization(miRNAdata = df,
                              report_summary = FALSE)

  # Whether the best_method is valid
  expect_true(res$best_method %in% names(res$normalized))
})


test_that("Check if normalized values are reasonable", {

  set.seed(624)
  df <- matrix(abs(rnorm(20, mean = 5)), nrow = 5, ncol = 4)
  res <- compareNormalization(miRNAdata = df,
                              report_summary = FALSE)

  # log2
  expect_true(all(res$normalized$log2 >= 0))

  # z-score
  expect_equal(unname(round(colMeans(res$normalized$zscore), 6)),
               rep(0, ncol(df)))

  # quantile normalization: columns should have same sorted values
  sorted_cols <- apply(res$normalized$quantile, 2, sort)
  expect_true(all(apply(sorted_cols, 1, function(x) length(unique(x)) == 1)))
})


test_that("Check if edge cases are correctly handled", {

  set.seed(624)
  df <- matrix(5, nrow = 5, ncol = 4)
  expect_warning(compareNormalization(miRNAdata = df,
                                      report_summary = FALSE),
                 "Some columns have zero variance")

  # miRNASeq1 also contains constant columns, warning should appear
  expect_warning(res <- compareNormalization(miRNAdata = miRNASeq1,
                                             report_summary = FALSE),
                 "Some columns have zero variance")

  # Whether the process still continues
  expect_true(all(c("normalized", "best_method") %in% names(res)))
  expect_true(all(c("log2", "zscore", "quantile") %in% names(res$normalized)))

  # Whether dataset dimentions stay the same
  expect_equal(dim(res$normalized$log2), dim(miRNASeq1))
})


test_that("Check if compareNormalization is deterministic", {

  set.seed(624)
  df <- abs(matrix(rnorm(20), nrow = 5, ncol = 4))
  res1 <- compareNormalization(miRNAdata = df, report_summary = FALSE)
  res2 <- compareNormalization(miRNAdata = df, report_summary = FALSE)

  # Whether results don't change when process repeated
  expect_equal(res1$normalized$zscore, res2$normalized$zscore)
  expect_equal(res1$best_method, res2$best_method)
})

