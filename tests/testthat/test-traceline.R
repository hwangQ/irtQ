# ---- helpers -------------------------------------------------------------------

# item metadata with a 3PLM item, a GRM item, and a GPCM item
make_mixed_meta <- function() {
  data.frame(
    id = c("D3", "G4", "P5"), cats = c(2, 4, 5), model = c("3PLM", "GRM", "GPCM"),
    par.1 = c(1.2, 1.4, 0.9), par.2 = c(0.4, -1.2, -1), par.3 = c(0.18, 0.1, 0),
    par.4 = c(NA, 1.3, 0.7), par.5 = c(NA, NA, 1.6)
  )
}

# ---- estimation objects --------------------------------------------------------

test_that("traceline() for est_irt and est_item objects uses the stored metadata and D", {
  meta <- make_mixed_meta()
  theta <- seq(-3, 3, 1)
  ref <- traceline(meta, theta, D = 1.702)
  fit_irt <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_irt", "list"))
  fit_item <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_item", "list"))
  expect_identical(traceline(fit_irt, theta), ref)
  expect_identical(traceline(fit_item, theta), ref)

  # a D supplied through ... does not replace the stored scaling constant
  expect_identical(traceline(fit_irt, theta, D = 1), ref)
})

# ---- structure of the results ---------------------------------------------------

test_that("traceline() keeps one-row matrices for a single theta", {
  meta <- make_mixed_meta()
  tr <- traceline(meta, 0)
  expect_true(all(vapply(tr$prob.cats, is.matrix, logical(1))))
  expect_equal(vapply(tr$prob.cats, ncol, numeric(1)), c(D3 = 2, G4 = 4, P5 = 5))
  expect_equal(vapply(tr$prob.cats, nrow, numeric(1)), c(D3 = 1, G4 = 1, P5 = 1))

  # the test characteristic curve can be plotted for a single theta
  grDevices::pdf(NULL)
  p <- plot(tr)
  grDevices::dev.off()
  expect_s3_class(p, "ggplot")
})

test_that("traceline() ICC and TCC are expected scores of the category probabilities", {
  meta <- make_mixed_meta()
  tr <- traceline(meta, seq(-3, 3, 1))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-8)
  expect_equal(tr$tcc, rowSums(tr$icc))
  expect_true(all(abs(sapply(tr$prob.cats, rowSums) - 1) < 1e-8))
})

test_that("traceline() ICC ignores the categories an item does not have", {
  # a GRM item with 3 categories next to a GRM item with 5 categories
  meta <- data.frame(
    id = c("A", "B"), cats = c(3, 5), model = "GRM", par.1 = c(1, 1.2), par.2 = c(-1, -1),
    par.3 = c(1, 0), par.4 = c(NA, 1), par.5 = c(NA, 2)
  )
  tr <- traceline(meta, c(-2, 0, 1.5))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-12)
  expect_equal(tr$tcc, rowSums(es), tolerance = 1e-12)

  # the same holds for GPCM items
  meta$model <- "GPCM"
  tr <- traceline(meta, c(-2, 0, 1.5))
  es <- sapply(tr$prob.cats, function(p) p %*% (0:(ncol(p) - 1)))
  expect_equal(unname(tr$icc), unname(es), tolerance = 1e-12)
})

test_that("traceline() computes a GRM or GPCM item with two categories as a 2PLM item", {
  theta <- c(-2, -1, 0, 0.5, 1.5)
  two_cat <- data.frame(
    id = c("G", "P"), cats = 2, model = c("GRM", "GPCM"), par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  tr <- traceline(two_cat, theta, D = 1.702)

  # the probabilities equal those of prm() with a single threshold
  expect_equal(tr$prob.cats$G[, 2], prm(theta, 1.2, 0.3, 1.702, "GRM")[, 2], tolerance = 1e-12)
  expect_equal(tr$prob.cats$P[, 2], prm(theta, 0.8, -0.4, 1.702, "GPCM")[, 2], tolerance = 1e-12)

  # and the 2PLM items with the same parameters
  twopl <- data.frame(
    id = c("G", "P"), cats = 2, model = "2PLM", par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  expect_equal(tr$icc, traceline(twopl, theta, D = 1.702)$icc, tolerance = 1e-12)
})

test_that("traceline() treats NA guessing parameters of 3PLM and DRM items as zeros", {
  na_g <- data.frame(
    id = c("A", "B"), cats = 2, model = c("3PLM", "DRM"), par.1 = 1, par.2 = 0, par.3 = NA
  )
  tr <- suppressWarnings(traceline(na_g, c(-1, 0, 1)))
  expect_equal(unname(tr$icc[, 1]), plogis(c(-1, 0, 1)), tolerance = 1e-12)
  expect_equal(unname(tr$icc[, 2]), plogis(c(-1, 0, 1)), tolerance = 1e-12)
})
