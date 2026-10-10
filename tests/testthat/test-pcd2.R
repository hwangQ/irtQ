# Independent implementation of the pseudo-count D2 (Lim & Han, 2026, Eqs. 2-5)
pcd2_ref <- function(par, resp, D, quad, w) {
  w <- w / sum(w)
  g <- ifelse(is.na(par$par.3), 0, par$par.3)
  pmat <- sapply(seq_len(nrow(par)), function(j) {
    g[j] + (1 - g[j]) / (1 + exp(-D * par$par.1[j] * (quad - par$par.2[j])))
  })
  post <- t(apply(resp, 1, function(u) {
    lik <- rep(1, length(quad))
    for (j in which(!is.na(u))) lik <- lik * if (u[j] == 1) pmat[, j] else 1 - pmat[, j]
    lik * w / sum(lik * w)
  }))
  sapply(seq_len(nrow(par)), function(j) {
    obs <- !is.na(resp[, j])
    r1 <- colSums(post[obs & resp[, j] == 1, , drop = FALSE])
    r0 <- colSums(post[obs & resp[, j] == 0, , drop = FALSE])
    sum((r0 + r1) * (r1 / (r0 + r1) - pmat[, j])^2) / sum(obs)
  })
}

pcd2_sim <- function(seed = 21) {
  set.seed(seed)
  nitem <- 8
  par <- shape_df(
    par.drm = list(a = runif(nitem, 0.8, 2), b = rnorm(nitem), g = runif(nitem, 0, 0.2)),
    cats = 2, model = "3PLM"
  )
  resp <- simdat(par, theta = rnorm(300), D = 1)
  resp[sample(length(resp), 0.2 * length(resp))] <- NA
  list(par = par, resp = resp)
}

test_that("pcd2() purification works down to the last item and stops when no item is left", {
  sim <- pcd2_sim()
  stat <- pcd2(sim$par[1:3, ], sim$resp[, 1:3])$no_purify$ipd_stat$pcd2
  # a critical value that flags the two largest values starts the purification with
  # one or two items left to analyze
  rst <- suppressWarnings(pcd2(sim$par[1:3, ], sim$resp[, 1:3], crit.val = sort(stat)[2],
                               purify = TRUE, max.iter = 5, verbose = FALSE))
  expect_false(anyNA(rst$with_purify$ipd_stat$pcd2))
  expect_gte(rst$with_purify$n.iter, 1)
  # every item is flagged, so no item is left to analyze
  expect_error(
    pcd2(sim$par[1:3, ], sim$resp[, 1:3], crit.val = 1e-10, purify = TRUE, verbose = FALSE),
    "no item is left"
  )
})

test_that("pcd2() warns about min.resp only when an item is excluded", {
  sim <- pcd2_sim()
  expect_no_warning(pcd2(sim$par, sim$resp, min.resp = 5))
  resp <- sim$resp
  resp[-(1:20), 2] <- NA
  expect_warning(rst <- pcd2(sim$par, resp, min.resp = 50), "Item\\(s\\) 2 are not analyzed")
  expect_true(is.nan(rst$no_purify$ipd_stat$pcd2[2]))
})

test_that("pcd2() reads character responses and stops for invalid ones", {
  sim <- pcd2_sim()
  resp_chr <- sim$resp
  storage.mode(resp_chr) <- "character"
  expect_identical(pcd2(sim$par, resp_chr)$no_purify$ipd_stat, pcd2(sim$par, sim$resp)$no_purify$ipd_stat)
  resp_bad <- sim$resp
  resp_bad[1, 1] <- -1
  expect_error(pcd2(sim$par, resp_bad), "score categories")
  resp_bad[1, 1] <- 0.5
  expect_error(pcd2(sim$par, resp_bad), "score categories")
  expect_error(pcd2(sim$par, sim$resp[, -1]), "number of columns")
  expect_error(pcd2(sim$par, sim$resp, item.skip = "V2"), "item.skip")
})

test_that("pcd2() reads a logical item.skip as item positions", {
  sim <- pcd2_sim()
  skip_lgl <- seq_len(nrow(sim$par)) %in% c(2, 5)
  r1 <- pcd2(sim$par, sim$resp, item.skip = c(2, 5))
  r2 <- pcd2(sim$par, sim$resp, item.skip = skip_lgl)
  expect_identical(r1$no_purify, r2$no_purify)
  expect_true(all(is.na(r1$no_purify$ipd_stat$pcd2[c(2, 5)])))
  expect_error(pcd2(sim$par, sim$resp, item.skip = 20), "item.skip")
})
