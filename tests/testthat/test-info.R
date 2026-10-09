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

# ---- dichotomous information ---------------------------------------------------

test_that("info() of dichotomous items stays finite far from the item location", {
  meta <- data.frame(
    id = c("A", "B"), cats = 2, model = c("2PLM", "3PLM"), par.1 = 3, par.2 = 0,
    par.3 = c(NA, 0.2)
  )
  iif <- info(meta, c(-300, 300, 1000), D = 1.702)$iif
  expect_false(anyNA(iif))
  expect_equal(unname(iif), matrix(0, 2, 3), tolerance = 1e-8)
})

test_that("info() of dichotomous items matches the closed-form 3PLM information", {
  meta <- data.frame(id = c("A", "B"), cats = 2, model = c("3PLM", "2PLM"),
                     par.1 = c(1.2, 0.8), par.2 = c(0.4, -0.5), par.3 = c(0.18, NA))
  theta <- seq(-3, 3, 1)
  for (D in c(1, 1.702)) {
    iif <- info(meta, theta, D = D)$iif
    # information (D a)^2 (Q / P) ((P - g) / (1 - g))^2 with P = g + (1 - g) / (1 + exp(-D a (theta - b)))
    ref <- function(a, b, g) {
      P <- g + (1 - g) / (1 + exp(-D * a * (theta - b)))
      (D * a)^2 * ((1 - P) / P) * ((P - g) / (1 - g))^2
    }
    expect_equal(unname(iif[1, ]), ref(1.2, 0.4, 0.18), tolerance = 1e-8)
    expect_equal(unname(iif[2, ]), ref(0.8, -0.5, 0), tolerance = 1e-8)
  }
})

test_that("info() computes a GRM or GPCM item with two categories as a 2PLM item", {
  theta <- c(-2, -1, 0, 0.5, 1.5)
  two_cat <- data.frame(
    id = c("G", "P"), cats = 2, model = c("GRM", "GPCM"), par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  twopl <- data.frame(
    id = c("G", "P"), cats = 2, model = "2PLM", par.1 = c(1.2, 0.8),
    par.2 = c(0.3, -0.4), par.3 = NA
  )
  expect_equal(
    info(two_cat, theta, D = 1.702)$iif,
    info(twopl, theta, D = 1.702)$iif,
    tolerance = 1e-12
  )
})

# ---- plot labels ---------------------------------------------------------------

test_that("plot.info() labels items with their IDs as given", {
  meta <- data.frame(
    id = c("1", "item-2", "CR 3"), cats = c(2, 2, 3), model = c("3PLM", "2PLM", "GRM"),
    par.1 = c(1, 1.2, 1), par.2 = c(0, 0.5, -1), par.3 = c(0.2, NA, 1)
  )
  inf <- info(meta, seq(-2, 2, 0.5))
  p <- plot(inf, item.loc = 1:3)
  expect_setequal(as.character(unique(p$data$item)), c("1", "item-2", "CR 3"))
  expect_identical(plot(inf, item.loc = 2)$labels$title, "Item Information: item-2")
})
