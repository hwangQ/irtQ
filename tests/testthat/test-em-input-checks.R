# Input checks and item IDs of est_irt()

xf <- bring.flexmirt(file = system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ"), "par")$Group1$full_df
x12 <- xf[c(1:10, 39:40), ]
set.seed(1)
d12 <- simdat(x = x12, theta = rnorm(500), D = 1)

test_that("est_irt() stops with a clear message on invalid responses", {
  dbad <- d12
  dbad[1, 1] <- 2
  msg <- paste0("scored as 0, 1.*", x12$id[1], " \\(allowed 0 to 1\\)")
  expect_error(est_irt(x = x12, data = dbad, D = 1, verbose = FALSE), msg)
  expect_error(est_irt(x = x12, data = dbad, D = 1, fipc = TRUE, fix.loc = 1:5, verbose = FALSE), msg)
  dneg <- d12
  dneg[2, 11] <- -1
  expect_error(est_irt(x = x12, data = dneg, D = 1, verbose = FALSE),
               paste0(x12$id[11], " \\(allowed 0 to 4\\)"))
  dfrac <- d12
  dfrac[3, 2] <- 0.5
  expect_error(est_irt(x = x12, data = dfrac, D = 1, verbose = FALSE), "non-integer")
})

test_that("est_irt() with FIPC stops with a clear message on invalid inputs", {
  fipc_fit <- function(...) est_irt(x = x12, D = 1, fipc = TRUE, verbose = FALSE, ...)
  expect_error(fipc_fit(data = cbind(d12, d12[, 1]), fix.loc = 1:5), "must be the same")
  expect_error(fipc_fit(data = d12[, -12], fix.loc = 1:5), "must be the same")
  expect_error(fipc_fit(data = d12, fix.loc = c(1, 20)), "'fix.loc' must contain row positions of 'x' between 1 and 12")
  expect_error(fipc_fit(data = d12, fix.loc = integer(0)), "'fix.loc' must contain row positions")
  expect_error(fipc_fit(data = d12, fix.loc = c(1:3, 1)), "duplicated")
  expect_error(fipc_fit(data = d12, fix.id = "none"), "'fix.id' not found.*none")
  expect_error(fipc_fit(data = d12, fix.id = c(x12$id[1:2], "none")), "'fix.id' not found.*none")
})

test_that("item.id replaces the IDs in x with and without FIPC", {
  new.id <- paste0("Z", seq_len(nrow(x12)))
  fit <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, item.id = new.id, MaxE = 3, se = FALSE, verbose = FALSE))
  expect_identical(fit$par.est$id, new.id)
  expect_identical(fit$estimates$id, new.id)
  fit.fipc <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, item.id = new.id, fipc = TRUE, fix.id = new.id[1:5],
                                       MaxE = 3, se = FALSE, verbose = FALSE))
  expect_identical(fit.fipc$par.est$id, new.id)
  expect_error(est_irt(x = x12, data = d12, D = 1, item.id = new.id[-1], verbose = FALSE), "length of 'item.id'")
  expect_error(est_irt(x = x12, data = d12, D = 1, item.id = new.id[-1], fipc = TRUE, fix.loc = 1:5, verbose = FALSE),
               "length of 'item.id'")
})

test_that("FIPC with all items fixed reports that no information matrix is computed", {
  fit <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, fipc = TRUE, fix.loc = seq_len(nrow(x12)), verbose = FALSE))
  expect_identical(fit$test.2, "No item parameters are estimated, so the information matrix is not computed.")
  # the message of a FIPC fit with items to estimate is unchanged
  fit2 <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, fipc = TRUE, fix.loc = 1:6, verbose = FALSE))
  expect_match(fit2$test.2, "possible local maximum")
})

test_that("FIPC with all items fixed compares the first cycle with the moments of the supplied weights", {
  w <- gen.weight(dist = "norm", mu = 0.5, sigma = 1.2, theta = seq(-6, 6, length.out = 49))
  fit <- suppressWarnings(est_irt(x = x12, data = d12, D = 1, fipc = TRUE, fix.loc = seq_len(nrow(x12)), weights = w,
                                  MaxE = 1, verbose = FALSE))
  mom.w <- cal_moment(node = w[, 1], weight = w[, 2])
  mom.new <- cal_moment(node = fit$weights[, 1], weight = fit$weights[, 2])
  expect_equal(fit$maxpar.diff, max(abs(mom.new - mom.w)), tolerance = 1e-12)
})
