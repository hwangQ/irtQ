# convergence of the EM algorithm is judged by the largest absolute change in all estimated parameters

# largest absolute change in the item parameters and in the listed latent distribution moments between
# two fits; 'mom' returns the mean and variance of the distributions that are estimated as a numeric vector
em_change <- function(f1, f0, items, mom) {
  d_item <- abs(as.matrix(items(f1)) - as.matrix(items(f0)))
  max(c(d_item, abs(mom(f1) - mom(f0))), na.rm = TRUE)
}

irt_items <- function(f) f$par.est[, grepl("^par\\.", names(f$par.est)), drop = FALSE]
irt_mom <- function(f) unlist(f$group.par[1, c("mu", "sigma2")])
mg_items <- function(f) f$par.est$overall[, grepl("^par\\.", names(f$par.est$overall)), drop = FALSE]
mg_mom <- function(free) {
  function(f) unlist(lapply(free, function(g) unlist(f$group.par[[g]]["estimates", c("mu", "sigma2")])))
}

test_that("est_irt with FIPC stops on the change in the items and in the prior mean and variance", {
  skip_on_cran()
  set.seed(21)
  x <- shape_df(par.drm = list(a = rlnorm(10, 0, .25), b = rnorm(10), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(400, 0.3, 1.2), D = 1)
  for (emp in c(FALSE, TRUE)) {
    # the histogram moments move slowly, so the empirical histogram uses a looser criterion
    etol <- if (emp) 1e-3 else 1e-4
    fit <- function(MaxE) {
      # a fit that is stopped before convergence warns about it
      suppressWarnings(est_irt(x = x, data = d, D = 1, fipc = TRUE, fix.loc = 1:5, EmpHist = emp, Etol = etol,
                               MaxE = MaxE, se = FALSE, verbose = FALSE))
    }
    f <- fit(500)
    expect_lte(f$maxpar.diff, etol)
    f0 <- fit(f$niter - 1)
    expect_equal(em_change(f, f0, irt_items, irt_mom), f$maxpar.diff, tolerance = 1e-8)
  }
})

test_that("est_irt without FIPC stops on the change in the item parameters", {
  skip_on_cran()
  set.seed(22)
  x <- shape_df(par.drm = list(a = rlnorm(8, 0, .25), b = rnorm(8), g = NA), cats = 2, model = "2PLM")
  d <- simdat(x = x, theta = rnorm(500), D = 1)
  fit <- function(MaxE) {
    suppressWarnings(est_irt(data = d, D = 1, model = "2PLM", cats = 2, MaxE = MaxE, se = FALSE, verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  expect_equal(em_change(f, f0, irt_items, function(f) numeric(0)), f$maxpar.diff, tolerance = 1e-8)
})

test_that("est_mg stops on the change in the items and in the means and variances of the free groups", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  for (emp in c(FALSE, TRUE)) {
    fit <- function(MaxE) {
      suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, EmpHist = emp,
                              MaxE = MaxE, se = FALSE, verbose = FALSE))
    }
    f <- fit(500)
    expect_lte(f$maxpar.diff, f$Etol)
    f0 <- fit(f$niter - 1)
    expect_equal(em_change(f, f0, mg_items, mg_mom(2)), f$maxpar.diff, tolerance = 1e-8)
  }
})

test_that("est_mg with FIPC stops on the change in the free items and in the free groups", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  fit <- function(MaxE) {
    suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, fipc = TRUE,
                            fix.id = c("C1", "C2"), MaxE = MaxE, se = FALSE, verbose = FALSE))
  }
  f <- fit(500)
  expect_lte(f$maxpar.diff, f$Etol)
  f0 <- fit(f$niter - 1)
  expect_equal(em_change(f, f0, mg_items, mg_mom(2)), f$maxpar.diff, tolerance = 1e-8)
})

test_that("a group that is not in free.group does not enter the convergence statistic", {
  skip_on_cran()
  M <- make_mg_data(seed = 11)
  gn <- c("G1", "G2")
  # all items fixed: only the distribution of the free group is estimated
  f <- suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, fipc = TRUE,
                               fix.id = M$xall$id[1:14], EmpHist = TRUE, se = FALSE, verbose = FALSE))
  f0 <- suppressWarnings(est_mg(x = M$x[1:2], data = M$data[1:2], group.name = gn, free.group = 2, fipc = TRUE,
                                fix.id = M$xall$id[1:14], EmpHist = TRUE, MaxE = f$niter - 1, se = FALSE, verbose = FALSE))
  expect_lte(f$maxpar.diff, f$Etol)
  expect_equal(em_change(f, f0, mg_items, mg_mom(2)), f$maxpar.diff, tolerance = 1e-8)
})
