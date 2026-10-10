# Independent implementation of the RIPD statistics (Lim & Han, 2026, Eqs. 6-17)
ripd_ref <- function(par, resp, score, group, focal, D) {
  out <- NULL
  g <- ifelse(is.na(par$par.3), 0, par$par.3)
  for (j in seq_len(nrow(par))) {
    P <- g[j] + (1 - g[j]) / (1 + exp(-D * par$par.1[j] * (score - par$par.2[j])))
    u <- resp[, j]
    ok <- !is.na(u) & !is.na(score)
    f <- ok & group == focal
    r <- ok & group != focal
    e <- u - P
    q <- P * (1 - P)
    vr <- mean(e[f]) - mean(e[r])
    vs <- mean(e[f]^2) - mean(e[r]^2)
    mus <- mean(q[f]) - mean(q[r])
    s2r <- sum(q[f]) / sum(f)^2 + sum(q[r]) / sum(r)^2
    s2s <- sum(q[f] * (1 - 2 * P[f])^2) / sum(f)^2 + sum(q[r] * (1 - 2 * P[r])^2) / sum(r)^2
    crs <- sum(q[f] * (1 - 2 * P[f])) / sum(f)^2 + sum(q[r] * (1 - 2 * P[r])) / sum(r)^2
    dv <- c(vr, vs - mus)
    chi <- drop(dv %*% solve(matrix(c(s2r, crs, crs, s2s), 2), dv))
    out <- rbind(out, data.frame(
      ripdr = vr, z.ripdr = vr / sqrt(s2r), ripds = vs, z.ripds = (vs - mus) / sqrt(s2s),
      ripdrs = chi, p.ripdr = 2 * pnorm(-abs(vr / sqrt(s2r))),
      p.ripds = 2 * pnorm(-abs((vs - mus) / sqrt(s2s))), p.ripdrs = pchisq(chi, 2, lower.tail = FALSE),
      n.ref = sum(r), n.foc = sum(f), mu.ripds = mus, sigma.ripdr = sqrt(s2r),
      sigma.ripds = sqrt(s2s), covariance = crs
    ))
  }
  out
}

# Simulate a small two-group data set with drift on the first three items
ripd_sim <- function(seed = 11) {
  set.seed(seed)
  nitem <- 12
  par <- shape_df(
    par.drm = list(a = runif(nitem, 0.8, 2), b = rnorm(nitem), g = c(rep(0, 4), runif(nitem - 4, 0, 0.25))),
    cats = 2, model = c(rep("2PLM", 4), rep("3PLM", nitem - 4))
  )
  par_drift <- par
  par_drift$par.2[1:3] <- par_drift$par.2[1:3] + 0.6
  resp <- rbind(simdat(par, theta = rnorm(600), D = 1), simdat(par_drift, theta = rnorm(400, 0.2), D = 1))
  resp[sample(length(resp), 0.25 * length(resp))] <- NA
  resp <- resp[rowSums(!is.na(resp)) > 0, ]
  group <- c(rep("r", 600), rep("f", 400))[seq_len(nrow(resp))]
  list(par = par, resp = resp, group = group)
}

test_that("ripd() excludes examinees with a missing ability estimate with a warning", {
  sim <- ripd_sim()
  score <- suppressWarnings(est_score(sim$par, sim$resp, D = 1)$est.theta)
  score[c(1:30, 701:720)] <- NA
  expect_warning(
    rst <- ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f", verbose = FALSE),
    "missing ability estimate"
  )
  ref <- ripd_ref(sim$par, sim$resp, score, sim$group, "f", 1)
  expect_equal(rst$no_purify$ipd_stat$ripdrs, ref$ripdrs, tolerance = 1e-4)
  expect_equal(rst$no_purify$ipd_stat$n.ref, ref$n.ref)
  expect_equal(rst$no_purify$ipd_stat$n.foc, ref$n.foc)
})

test_that("ripd() applies min.resp to the focal group before scoring", {
  sim <- ripd_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- sim$group == "f" & n_resp < 8
  resp <- sim$resp
  resp[drop_row, ] <- NA
  r1 <- suppressWarnings(ripd(sim$par, sim$resp, group = sim$group, focal.name = "f", min.resp = 8, verbose = FALSE))
  r2 <- suppressWarnings(ripd(sim$par, resp, group = sim$group, focal.name = "f", verbose = FALSE))
  expect_true(any(drop_row))
  expect_identical(r1$no_purify$ipd_stat, r2$no_purify$ipd_stat)
})

test_that("ripd() applies min.resp to the focal group when score is provided", {
  sim <- ripd_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- sim$group == "f" & n_resp < 8
  score <- suppressWarnings(est_score(sim$par, sim$resp, D = 1)$est.theta)
  score_na <- score
  score_na[drop_row] <- NA
  resp_na <- sim$resp
  resp_na[drop_row, ] <- NA
  expect_true(any(drop_row))
  # the examinees below min.resp are excluded with a warning about their number
  expect_warning(
    r1 <- ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f",
               min.resp = 8, verbose = FALSE),
    paste(sum(drop_row), "focal group examinee")
  )
  r2 <- suppressWarnings(ripd(sim$par, resp_na, score = score_na, group = sim$group,
                              focal.name = "f", verbose = FALSE))
  expect_identical(r1$no_purify$ipd_stat, r2$no_purify$ipd_stat)
})

test_that("ripd() purification works when items are removed down to the last one", {
  x1 <- shape_df(par.drm = list(a = c(1, 1.2, 0.9), b = c(0, 0.5, -0.5), g = 0), cats = 2, model = "2PLM")
  x2 <- x1
  x2$par.2 <- x2$par.2 + c(1.5, 1.2, 0)
  set.seed(5)
  resp <- rbind(simdat(x1, theta = rnorm(500), D = 1), simdat(x2, theta = rnorm(500), D = 1))
  group <- rep(0:1, each = 500)
  rst <- ripd(x1, resp, group = group, focal.name = 1, purify = TRUE, verbose = FALSE)
  expect_true(length(rst$with_purify$ipd_item) >= 2)
  expect_false(anyNA(rst$with_purify$ipd_stat$ripdrs))
  rst2 <- ripd(x1[1:2, ], resp[, 1:2], group = group, focal.name = 1, purify = TRUE, verbose = FALSE)
  expect_true(rst2$with_purify$complete)
})

test_that("ripd() reads a logical item.skip as item positions during purification", {
  sim <- ripd_sim()
  skip_num <- c(1, 11)
  skip_lgl <- seq_len(nrow(sim$par)) %in% skip_num
  r1 <- suppressWarnings(ripd(sim$par, sim$resp, group = sim$group, focal.name = "f",
                              item.skip = skip_num, purify = TRUE, verbose = FALSE))
  r2 <- suppressWarnings(ripd(sim$par, sim$resp, group = sim$group, focal.name = "f",
                              item.skip = skip_lgl, purify = TRUE, verbose = FALSE))
  expect_identical(r1$with_purify, r2$with_purify)
  expect_true(all(is.na(r1$with_purify$ipd_stat[skip_num, 2:9])))
})

test_that("ripd() stops when purification flags every item", {
  # alpha = 1 flags every item at every iteration
  x1 <- shape_df(par.drm = list(a = c(1, 1.2), b = c(0, 0.5), g = 0), cats = 2, model = "2PLM")
  x2 <- x1
  x2$par.2 <- x2$par.2 + c(2, -2)
  set.seed(5)
  resp <- rbind(simdat(x1, theta = rnorm(500), D = 1), simdat(x2, theta = rnorm(500), D = 1))
  group <- rep(0:1, each = 500)
  expect_error(
    ripd(x1, resp, group = group, focal.name = 1, purify = TRUE, alpha = 1, verbose = FALSE),
    "no item is left"
  )
})

test_that("ripd() reads character responses as numbers", {
  sim <- ripd_sim()
  score <- suppressWarnings(est_score(sim$par, sim$resp, D = 1)$est.theta)
  resp_chr <- sim$resp
  storage.mode(resp_chr) <- "character"
  r1 <- ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f", verbose = FALSE)
  r2 <- ripd(sim$par, resp_chr, score = score, group = sim$group, focal.name = "f", verbose = FALSE)
  r3 <- ripd(sim$par, as.data.frame(resp_chr), score = score, group = sim$group, focal.name = "f", verbose = FALSE)
  expect_identical(r1$no_purify$ipd_stat, r2$no_purify$ipd_stat)
  expect_identical(r1$no_purify$ipd_stat, r3$no_purify$ipd_stat)
})

test_that("ripd() stops with clear errors for invalid inputs", {
  sim <- ripd_sim()
  score <- suppressWarnings(est_score(sim$par, sim$resp, D = 1)$est.theta)
  n <- nrow(sim$resp)
  resp_bad <- sim$resp
  resp_bad[1, 1] <- -1
  expect_error(ripd(sim$par, resp_bad, score = score, group = sim$group, focal.name = "f"), "score categories")
  resp_bad[1, 1] <- 0.5
  expect_error(ripd(sim$par, resp_bad, score = score, group = sim$group, focal.name = "f"), "score categories")
  resp_bad <- sim$resp
  resp_bad[is.na(resp_bad)] <- -9
  expect_error(ripd(sim$par, resp_bad, score = score, group = sim$group, focal.name = "f"), "score categories")
  expect_error(ripd(sim$par, sim$resp[, -1], score = score, group = sim$group, focal.name = "f"), "number of columns")
  expect_error(ripd(sim$par, sim$resp, score = score[-1], group = sim$group, focal.name = "f"), "'score'")
  expect_error(ripd(sim$par, sim$resp, score = score, group = sim$group[-1], focal.name = "f"), "'group'")
  expect_error(ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "F"), "focal.name")
  expect_error(ripd(sim$par, sim$resp, score = score, group = rep("f", n), focal.name = "f"), "reference group")
  expect_error(ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = 20), "item.skip")
  expect_error(ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f", item.skip = "V2"), "item.skip")
  xp <- shape_df(par.prm = list(a = 1, d = list(c(0, 1))), cats = 3, model = "GPCM")
  expect_error(ripd(rbind(sim$par, xp), cbind(sim$resp, 0), score = score, group = sim$group, focal.name = "f"),
               "dichotomous")
})

test_that("ripd() flags items with the unrounded p-values", {
  sim <- ripd_sim()
  score <- suppressWarnings(est_score(sim$par, sim$resp, D = 1)$est.theta)
  ref <- ripd_ref(sim$par, sim$resp, score, sim$group, "f", 1)

  # choose an item whose p-value changes by rounding and put alpha between the two values
  p_round <- round(ref$p.ripdr, 4)
  j <- which(abs(p_round - ref$p.ripdr) > 1e-5 & ref$p.ripdr > 1e-3)[1]
  expect_false(is.na(j))
  alpha <- (p_round[j] + ref$p.ripdr[j]) / 2
  rst <- ripd(sim$par, sim$resp, score = score, group = sim$group, focal.name = "f", alpha = alpha,
              verbose = FALSE)
  expect_equal(j %in% rst$no_purify$ipd_item$ripdr, ref$p.ripdr[j] <= alpha)
  expect_equal(rst$no_purify$ipd_stat$p.ripdr[j], p_round[j])
})
