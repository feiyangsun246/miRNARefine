library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if miRNAStability error upon invalid dataset input", {

  # Handling missing values to avoid unexpected errors
  filled1 <- missingValueHandling(miRNAdata = miRNASeq1, method = "median",
                                  report_summary = FALSE)
  filled2 <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
                                  report_summary = FALSE)

  # Valid input
  expect_silent(miRNAStability(miRNAdata = filled2, report_summary = FALSE))

  # Invalid input: Vector
  expect_error(miRNAStability(miRNAdata = miRNASeq1[, 1]),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: List
  expect_error(miRNAStability(miRNAdata = as.list(miRNASeq1)),
               "`miRNAdata` must be a data frame or matrix")

  # Invalid input: Numeric
  expect_error(miRNAStability(miRNAdata = 624),
               "`miRNAdata` must be a data frame or matrix")

  # Empty input
  df_empty <- data.frame()
  expect_error(miRNAStability(miRNAdata = df_empty),
               "Empty dataframe input")
})


test_that("Check if output structure is correct", {

  set.seed(624)
  df <- matrix(runif(20, 1, 10), nrow = 5, ncol = 4)
  res <- miRNAStability(miRNAdata = df, report_summary = FALSE)

  # Whether result is a list
  expect_type(res, "list")

  # Whether result contains certain domains
  expect_true(all(c("stability_scores", "most_stable",
                    "least_stable") %in% names(res)))
  expect_s3_class(res$stability_scores, "data.frame")
  expect_true(all(c("miRNA", "CV", "MAD") %in% colnames(res$stability_scores)))
})


test_that("Check if stability metrics make sense", {

  # Constant columns
  df <- matrix(rep(1:5, each = 4), nrow = 5)
  res <- miRNAStability(miRNAdata = df, report_summary = FALSE)

  expect_true(all(res$stability_scores$CV >= 0))
  expect_true(all(res$stability_scores$MAD >= 0))
})


test_that("Check if correctly computes CV and MAD for simple data", {
  df <- matrix(rep(1:4, each = 3), nrow = 3)
  colnames(df) <- paste0("miR", 1:4)

  res <- miRNAStability(miRNAdata = df, report_summary = FALSE)

  # Compute manually
  expected_cv <- apply(df, 2, function(x) sd(x) / mean(x))

  expect_equal(round(res$stability_scores$CV, 6),
               round(as.numeric(expected_cv), 6))
})


test_that("Check if handles constant columns correctly", {

  # Elements all the same
  df <- matrix(5, nrow = 4, ncol = 3)
  res <- miRNAStability(miRNAdata = df, report_summary = FALSE)

  expect_true(all(is.na(res$stability_scores$CV) | res$stability_scores$CV == 0))
  expect_true(all(res$stability_scores$MAD == 0))
})


test_that("Check if handles NA values correctly", {

  set.seed(624)
  df <- matrix(runif(20), nrow = 5)
  df[1, 1] <- NA
  expect_error(miRNAStability(miRNAdata = df, report_summary = FALSE),
               "Dataset contains missing values")
})


test_that("Check if top stable miRNAs selected correctly on simple data", {

  set.seed(624)
  df <- matrix(runif(50, 1, 10), nrow = 10)
  res1 <- miRNAStability(miRNAdata = df, num_top = 3, report_summary = FALSE)
  res2 <- miRNAStability(miRNAdata = df, report_summary = FALSE)

  # Check res1
  expect_length(res1$most_stable, 3)
  expect_true(all(res1$most_stable %in% res1$stability_scores$miRNA))

  # Check res2
  expect_length(res2$most_stable, 5)
  expect_true(all(res2$most_stable %in% res1$stability_scores$miRNA))
})


test_that("Check if top stable miRNAs selected correctly on actual
          miRNA data", {

  res1 <- miRNAStability(miRNAdata = miRNASeq1, num_top = 3,
                         report_summary = FALSE)
  res2 <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)

  # Check res1
  expect_length(res1$most_stable, 3)
  expect_true(all(res1$most_stable %in% res1$stability_scores$miRNA))

  # Check res2
  expect_length(res2$most_stable, 5)
  expect_true(all(res2$most_stable %in% res1$stability_scores$miRNA))
})


test_that("Check if report_summary runs silently when FALSE", {

  set.seed(624)
  df <- matrix(runif(20, 1, 10), nrow = 5)
  expect_silent(miRNAStability(miRNAdata = df, report_summary = FALSE))
})

