library(miRNARefine)
data("miRNASeq1")
data("miRNASeq2")

test_that("Check if handles empty or invalid input correctly", {

  empty_lst <- list()
  empty_df <- data.frame()

  # When input is NULL
  expect_error(plotStabilityDistribution(NULL), "Input is NULL")
  # When input dataset is empty
  expect_error(plotStabilityDistribution(stability_results = empty_lst),
               "Empty list")
  expect_error(plotStabilityDistribution(stability_results = empty_df),
               "Empty list")
})


test_that("Check if handles invalid metric input correctly", {

  res <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)

  # Invalid metric: something other than the two options
  expect_error(plotStabilityDistribution(stability_results = res,
                                         metric = "WHAT"),
               "`metric` must be either 'CV' or 'MAD'")

  # Invalid metric: empty
  expect_error(plotStabilityDistribution(stability_results = res,
                                         metric = ""),
               "`metric` must be either 'CV' or 'MAD'")

  expect_silent(plotStabilityDistribution(stability_results = res))
})


test_that("Check if returns a ggplot object", {

  # Using miRNASeq1
  res <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)
  p <- plotStabilityDistribution(stability_results = res)

  # Whether the output is a ggplot
  expect_s3_class(p, "ggplot")

  # Using miRNASeq2
  # Perform missingValueHandling first to avoid errors
  filled <- missingValueHandling(miRNAdata = miRNASeq2, method = "median",
                                              report_summary = FALSE)
  res2 <- miRNAStability(miRNAdata = filled, report_summary = FALSE)
  p2 <- plotStabilityDistribution(stability_results = res2)

  # Whether the output is a ggplot
  expect_s3_class(p2, "ggplot")

})


test_that("Check if handles input with wrong structure", {

  set.seed(624)
  df <- data.frame(rnorm(20), nrow = 5, ncol = 4)

  # Whether the process stops
  expect_error(plotStabilityDistribution(df), "Wrong format for input list")
})


test_that("Check if plot contains correct layers", {

  res <- miRNAStability(miRNAdata = miRNASeq1, report_summary = FALSE)
  p <- plotStabilityDistribution(stability_results = res)

  # At least one
  layer_classes <- sapply(p$layers, function(x) class(x$geom)[1])
  expect_true(any(layer_classes == "GeomPoint"))
})

