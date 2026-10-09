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

test_that("info() for est_irt and est_item objects uses the stored metadata and D", {
  meta <- make_mixed_meta()
  theta <- seq(-3, 3, 1)
  ref <- info(meta, theta, D = 1.702)
  fit_irt <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_irt", "list"))
  fit_item <- structure(list(par.est = meta, scale.D = 1.702), class = c("est_item", "list"))
  expect_identical(info(fit_irt, theta), ref)
  expect_identical(info(fit_item, theta), ref)

  # a D supplied through ... does not replace the stored scaling constant
  expect_identical(info(fit_irt, theta, D = 1), ref)

  # tif = FALSE leaves the tif component empty
  expect_null(info(fit_irt, theta, tif = FALSE)$tif)
})

# ---- GPCM information ----------------------------------------------------------

test_that("info() of a GPCM item does not depend on the other theta values", {
  meta <- data.frame(
    id = "G", cats = 5, model = "GPCM", par.1 = 2, par.2 = -1, par.3 = 0,
    par.4 = 1, par.5 = 2
  )
  alone <- info(meta, 0, D = 1.702)$iif[1, 1]
  together <- info(meta, c(0, 50), D = 1.702)$iif[1, 1]
  expect_equal(together, alone, tolerance = 1e-12)
  expect_equal(alone, 3.643564, tolerance = 1e-6)
})

test_that("info() of a GPCM item is zero far from the item location", {
  meta <- data.frame(
    id = "G", cats = 5, model = "GPCM", par.1 = 2, par.2 = -1, par.3 = 0,
    par.4 = 1, par.5 = 2
  )
  iif <- info(meta, c(-300, 300), D = 1.702)$iif
  expect_false(anyNA(iif))
  expect_equal(unname(iif[1, ]), c(0, 0), tolerance = 1e-8)
})
