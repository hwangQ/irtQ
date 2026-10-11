# Independent implementation of the RDIF-CR statistics from first principles.
# For an item with K categories the residual vector is the one-hot response
# minus the category probabilities. Conditional on the ability estimates, its
# covariance is diag(p) - p p', the covariance of the squared residuals is
# (diag(p) - p p') * w w' with w = 1 - 2p, and the cross covariance of a raw
# residual j and a squared residual k is (diag(p) - p p')[j, k] * w[k].
# Binary items use only the residual of the correct answer category.
crt_probs <- function(x, theta, D) {
  lapply(seq_len(nrow(x)), function(j) {
    K <- x$cats[j]
    a <- x$par.1[j]
    if (K == 2) {
      g <- if (x$model[j] == "3PLM") x$par.3[j] else 0
      P <- g + (1 - g) / (1 + exp(-D * a * (theta - x$par.2[j])))
      return(cbind(1 - P, P))
    }
    d <- unlist(x[j, paste0("par.", 2:K)])
    if (x$model[j] == "GRM") {
      ps <- cbind(1, sapply(d, function(b) 1 / (1 + exp(-D * a * (theta - b)))), 0)
      return(ps[, 1:K] - ps[, 2:(K + 1)])
    }
    z <- cbind(0, t(apply(sapply(d, function(b) D * a * (theta - b)), 1, cumsum)))
    e <- exp(z - apply(z, 1, max))
    e / rowSums(e)
  })
}

# quadratic form with the Moore-Penrose inverse and the rank of the matrix
crt_quad <- function(S, d) {
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * 1e-10
  list(chi = sum(crossprod(ev$vectors[, keep, drop = FALSE], d)^2 / ev$values[keep]), rank = sum(keep))
}

crt_ref <- function(x, resp, score, group, focal, D) {
  P <- crt_probs(x, ifelse(is.na(score), 0, score), D)
  out <- NULL
  mom <- list()
  for (j in seq_len(nrow(x))) {
    K <- x$cats[j]
    idx <- if (K == 2) 2 else seq_len(K)
    ok <- !is.na(resp[, j]) & !is.na(score)
    parts <- lapply(list(f = ok & group == focal, r = ok & group != focal), function(s) {
      p <- P[[j]][s, idx, drop = FALSE]
      y <- outer(resp[s, j], 0:(K - 1), "==")[, idx, drop = FALSE] * 1
      e <- y - p
      w <- 1 - 2 * p
      n <- sum(s)
      list(R = colMeans(e), S = colMeans(e^2), mu = colMeans(p * (1 - p)), n = n,
           Crr = (diag(colSums(p), length(idx)) - crossprod(p)) / n^2,
           Css = (diag(colSums(p * w^2), length(idx)) - crossprod(p * w)) / n^2,
           Crs = (diag(colSums(p * w), length(idx)) - crossprod(p, p * w)) / n^2)
    })
    f <- parts$f
    r <- parts$r
    R <- f$R - r$R
    S <- f$S - r$S - (f$mu - r$mu)
    Crr <- f$Crr + r$Crr
    Css <- f$Css + r$Css
    Crs <- f$Crs + r$Crs
    Call <- rbind(cbind(Crr, Crs), cbind(t(Crs), Css))
    qr <- crt_quad(Crr, R)
    qs <- crt_quad(Css, S)
    qrs <- crt_quad(Call, c(R, S))
    dfr <- if (K == 2) 1 else K - 1
    dfs <- if (K == 2) 1 else K
    out <- rbind(out, data.frame(
      crdifr = qr$chi, crdifs = qs$chi, crdifrs = qrs$chi,
      rank.r = qr$rank, rank.s = qs$rank, rank.rs = qrs$rank,
      df.r = dfr, df.s = dfs, df.rs = dfr + dfs,
      p.crdifr = pchisq(qr$chi, dfr, lower.tail = FALSE),
      p.crdifs = pchisq(qs$chi, dfs, lower.tail = FALSE),
      p.crdifrs = pchisq(qrs$chi, dfr + dfs, lower.tail = FALSE),
      n.ref = r$n, n.foc = f$n
    ))
    mom[[j]] <- list(mu.s = f$mu - r$mu, Crr = Crr, Css = Css, Call = Call)
  }
  attr(out, "moments") <- mom
  out
}

# Simulate a mixed-format two-group data set with global DIF on items 4 and 6
crt_sim <- function(seed = 31, miss = 0.25) {
  set.seed(seed)
  x <- shape_df(
    par.drm = list(a = c(1.2, 0.8, 1.6), b = c(-0.5, 0, 0.7), g = c(NA, 0.2, NA)),
    par.prm = list(a = c(1.2, 0.9, 1.5, 1),
                   d = list(c(-1, 0, 1), c(-0.8, 0.2, 1.1, 1.6), c(-0.5, 0.7), c(-1.2, -0.1, 0.9))),
    item.id = paste0("C", 1:7), cats = c(2, 2, 2, 4, 5, 3, 4),
    model = c("2PLM", "3PLM", "2PLM", "GRM", "GPCM", "GRM", "GPCM")
  )
  xf <- x
  xf[4, c("par.1", "par.2", "par.4")] <- list(2, -0.3, 0.4)
  xf[6, c("par.2", "par.3")] <- list(-0.1, 0.2)
  resp <- rbind(simdat(x, theta = rnorm(700, 0, 1.2), D = 1), simdat(xf, theta = rnorm(600, 0, 1.2), D = 1))
  resp[sample(length(resp), miss * length(resp))] <- NA
  keep <- rowSums(!is.na(resp)) > 0
  group <- c(rep(0, 700), rep(1, 600))
  list(x = x, resp = resp[keep, ], group = group[keep])
}

# Evaluate an expression and return its value together with all warning messages
crt_catch <- function(expr) {
  msgs <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msgs <<- c(msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = msgs)
}

test_that("crdif() statistics match a first-principles implementation", {
  sim <- crt_sim()
  for (D in c(1, 1.702)) {
    score <- suppressWarnings(est_score(sim$x, sim$resp, D = D, method = "ML")$est.theta)
    rst <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1, D = D)
    ref <- crt_ref(sim$x, sim$resp, score, sim$group, 1, D)
    stat <- rst$no_purify$dif_stat
    for (v in c("crdifr", "crdifs", "crdifrs", "p.crdifr", "p.crdifs", "p.crdifrs")) {
      expect_lt(max(abs(stat[[v]] - ref[[v]])), 1e-4)
    }
    expect_equal(stat$n.ref, ref$n.ref)
    expect_equal(stat$n.foc, ref$n.foc)

    # the degrees of freedom are the ranks of the covariance matrices:
    # K - 1, K, and 2K - 1 for a polytomous item and 1, 1, and 2 for a binary item
    expect_equal(stat$df.crdifr, ifelse(sim$x$cats == 2, 1, sim$x$cats - 1))
    expect_equal(stat$df.crdifs, ifelse(sim$x$cats == 2, 1, sim$x$cats))
    expect_equal(stat$df.crdifrs, stat$df.crdifr + stat$df.crdifs)
    expect_equal(stat$df.crdifr, ref$df.r)
    expect_equal(stat$df.crdifs, ref$df.s)
    expect_equal(stat$df.crdifrs, ref$df.rs)
    mom <- attr(ref, "moments")
    for (j in seq_len(nrow(sim$x))) {
      expect_equal(unname(rst$no_purify$moments$mu.crdifs[[j]]), unname(mom[[j]]$mu.s), tolerance = 1e-8)
      expect_equal(unname(rst$no_purify$moments$cov.crdifr[[j]]), unname(mom[[j]]$Crr), tolerance = 1e-8)
      expect_equal(unname(rst$no_purify$moments$cov.crdifs[[j]]), unname(mom[[j]]$Css), tolerance = 1e-8)
      expect_equal(unname(rst$no_purify$moments$cov.crdifrs[[j]]), unname(mom[[j]]$Call), tolerance = 1e-8)
    }
    expect_equal(rst$no_purify$dif_item$crdifrs, which(ref$p.crdifrs <= 0.05))
  }
})

test_that("crdif() covariance of the raw categorical residuals has rank K - 1", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1)
  ref <- crt_ref(sim$x, sim$resp, score, sim$group, 1, 1)
  poly <- sim$x$cats > 2

  # the residuals of the score categories sum to zero for each examinee
  for (j in which(poly)) {
    cv <- rst$no_purify$moments$cov.crdifr[[j]]
    expect_lt(max(abs(cv %*% rep(1, ncol(cv)))), 1e-12)
  }
  expect_equal(ref$rank.r[poly], sim$x$cats[poly] - 1)
  expect_equal(ref$rank.s[poly], sim$x$cats[poly])
  expect_equal(ref$rank.rs[poly], 2 * sim$x$cats[poly] - 1)
})

test_that("crdif() and rdif() agree on all three statistics for binary items", {
  sim <- crt_sim()
  b <- sim$x$cats == 2
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  r1 <- rdif(sim$x[b, ], sim$resp[, b], score = score, group = sim$group, focal.name = 1)
  r2 <- crdif(sim$x[b, ], sim$resp[, b], score = score, group = sim$group, focal.name = 1)
  s1 <- r1$no_purify$dif_stat
  s2 <- r2$no_purify$dif_stat
  expect_equal(s2$crdifr, s1$z.rdifr^2, tolerance = 1e-3)
  expect_equal(s2$crdifs, s1$z.rdifs^2, tolerance = 1e-3)
  expect_equal(s2$crdifrs, s1$rdifrs, tolerance = 1e-6)
  expect_equal(s2$p.crdifr, s1$p.rdifr, tolerance = 1e-3)
  expect_equal(s2$p.crdifrs, s1$p.rdifrs, tolerance = 1e-6)
  expect_equal(s2$df.crdifr, rep(1, sum(b)))
  expect_equal(s2$df.crdifrs, rep(2, sum(b)))
})

test_that("crdif() uses a generalized inverse when a covariance matrix is singular", {
  sim <- crt_sim(miss = 0)
  b <- sim$x$cats == 2
  expect_warning(
    rst <- crdif(sim$x[b, ], sim$resp[, b], score = rep(1, nrow(sim$resp)), group = sim$group,
                 focal.name = 1),
    "generalized inverse"
  )
  stat <- rst$no_purify$dif_stat

  # with a common ability the squared residual is a linear function of the raw residual
  expect_equal(stat$crdifrs, stat$crdifr, tolerance = 1e-6)
  expect_equal(stat$df.crdifrs, rep(1, sum(b)))
  expect_true(all(stat$crdifrs >= 0))
})

test_that("crdif() keeps the nominal level of the tests under no DIF", {
  skip_on_cran()
  set.seed(2026)
  cats <- c(3, 4, 5, 4, 3)
  x <- shape_df(par.prm = list(a = runif(5, 0.8, 1.8), d = lapply(cats, function(k) sort(rnorm(k - 1)))),
                cats = cats, model = c("GRM", "GPCM", "GRM", "GPCM", "GRM"))
  group <- rep(0:1, each = 1000)
  rate <- replicate(150, {
    theta <- rnorm(2000)
    resp <- simdat(x, theta = theta, D = 1)
    st <- crdif(x, resp, score = theta, group = group, focal.name = 1)$no_purify$dif_stat
    c(r = mean(st$p.crdifr <= 0.05), rs = mean(st$p.crdifrs <= 0.05))
  })

  # the rates with the degrees of freedom K and 2K would be about 0.02 and 0.03
  expect_gt(mean(rate["r", ]), 0.03)
  expect_lt(mean(rate["r", ]), 0.075)
  expect_gt(mean(rate["rs", ]), 0.03)
  expect_lt(mean(rate["rs", ]), 0.075)
})

test_that("crdif() purification removes the item with the largest statistic when the p-values underflow", {
  skip_on_cran()
  set.seed(7)
  x <- shape_df(par.prm = list(a = rep(1.2, 6), d = lapply(1:6, function(i) c(-1, 0, 1))),
                cats = 4, model = "GRM")
  xf <- x
  xf[1, c("par.2", "par.3", "par.4")] <- xf[1, c("par.2", "par.3", "par.4")] + 2.4
  xf[2, c("par.2", "par.3", "par.4")] <- xf[2, c("par.2", "par.3", "par.4")] + 3.2
  resp <- rbind(simdat(x, theta = rnorm(5000), D = 1), simdat(xf, theta = rnorm(5000), D = 1))
  group <- rep(0:1, each = 5000)
  rst <- crdif(x, resp, group = group, focal.name = 1, purify = TRUE, verbose = FALSE)
  stat <- rst$no_purify$dif_stat

  # the p-values of the two items are zero in double precision, and item 2 has the larger statistic
  expect_equal(pchisq(stat$crdifrs[1:2], stat$df.crdifrs[1:2], lower.tail = FALSE), c(0, 0))
  expect_gt(stat$crdifrs[2], stat$crdifrs[1])

  # item 2 is removed first
  expect_equal(rst$with_purify$dif_stat$n.iter[2], 0)
  expect_equal(rst$with_purify$dif_stat$n.iter[1], 1)
})

test_that("crdif() methods for est_irt and est_item objects match the default method", {
  sim <- crt_sim(miss = 0)
  fit <- est_irt(data = sim$resp, D = 1, model = c("2PLM", "3PLM", "2PLM", "GRM", "GPCM", "GRM", "GPCM"),
                 cats = sim$x$cats, verbose = FALSE)
  r1 <- crdif(fit, group = sim$group, focal.name = 1)
  r2 <- crdif(fit$par.est, data = fit$data, group = sim$group, focal.name = 1, D = fit$scale.D)
  expect_identical(r1$no_purify, r2$no_purify)
  score <- est_score(sim$x, sim$resp, D = 1)$est.theta
  fit2 <- est_item(x = sim$x, data = sim$resp, score = score, D = 1, verbose = FALSE)
  r3 <- crdif(fit2, group = sim$group, focal.name = 1)
  r4 <- crdif(fit2$par.est, data = fit2$data, score = fit2$score, group = sim$group,
              focal.name = 1, D = fit2$scale.D)
  expect_identical(r3$no_purify, r4$no_purify)

  # the data and the scaling factor are taken from the object
  expect_error(crdif(fit, data = sim$resp, group = sim$group, focal.name = 1), "taken from the object")
  expect_error(crdif(fit2, D = 1.7, group = sim$group, focal.name = 1), "taken from the object")
})

test_that("crdif() leaves skipped items out of the tests and the purification", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  r0 <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1)
  r1 <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1, item.skip = 4)
  expect_true(all(is.na(r1$no_purify$dif_stat[4, 2:10])))
  expect_equal(r1$no_purify$dif_stat[-4, ], r0$no_purify$dif_stat[-4, ])
  expect_false(4 %in% unlist(r1$no_purify$dif_item))

  # a logical item.skip gives the same purification as the item positions
  p1 <- suppressWarnings(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1,
                               item.skip = 4, purify = TRUE, verbose = FALSE))
  p2 <- suppressWarnings(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1,
                               item.skip = seq_len(7) == 4, purify = TRUE, verbose = FALSE))
  expect_identical(p1$with_purify, p2$with_purify)
  expect_false(4 %in% p1$with_purify$dif_item)
})

test_that("crdif() excludes examinees with a missing ability estimate with a warning", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  score[c(1:20, 801:810)] <- NA
  expect_warning(
    rst <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1),
    "30 examinee\\(s\\) with a missing ability estimate"
  )
  ref <- crt_ref(sim$x, sim$resp, score, sim$group, 1, 1)
  expect_lt(max(abs(rst$no_purify$dif_stat$crdifrs - ref$crdifrs)), 1e-3)
  expect_equal(rst$no_purify$dif_stat$n.ref, ref$n.ref)
  expect_equal(rst$no_purify$dif_stat$n.foc, ref$n.foc)
  expect_error(
    crdif(sim$x, sim$resp, score = rep(NA_real_, nrow(sim$resp)), group = sim$group, focal.name = 1),
    "No examinee has an ability estimate"
  )

  # the excluded examinees stay excluded during purification
  out <- crt_catch(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1,
                         purify = TRUE, verbose = FALSE))
  expect_true(any(grepl("30 examinee\\(s\\) with a missing ability estimate", out$warnings)))
  expect_false(anyNA(out$value$with_purify$dif_stat$crdifrs))
})

test_that("crdif() reports the examinees who lose their ability estimates in the purification once", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)

  # five examinees respond only to the fourth item, which is flagged and removed in the purification
  resp <- sim$resp
  resp[5:9, ] <- NA
  resp[5:9, 4] <- 1
  out <- crt_catch(crdif(sim$x, resp, score = score, group = sim$group, focal.name = 1,
                         purify = TRUE, verbose = FALSE))
  expect_true(4 %in% out$value$with_purify$dif_item)
  msg <- grep("examinee\\(s\\) were excluded during the purification", out$warnings, value = TRUE)
  expect_length(msg, 1)
  expect_gte(as.numeric(sub(" .*", "", msg)), 5)
  expect_false(any(grepl("NA values are returned", out$warnings)))
})

test_that("crdif() purification removes the item with the smallest p-value first", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- suppressWarnings(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1,
                                purify = TRUE, verbose = FALSE))
  ref <- crt_ref(sim$x, sim$resp, score, sim$group, 1, 1)
  first <- which.min(ref$p.crdifrs)
  expect_true(first %in% rst$with_purify$dif_item)
  expect_equal(rst$with_purify$dif_stat$n.iter[first], 0)
  expect_true(rst$with_purify$complete)

  # the final scores are the estimates without the flagged items
  keep <- setdiff(seq_len(nrow(sim$x)), rst$with_purify$dif_item)
  sc <- suppressWarnings(est_score(sim$x[keep, ], sim$resp[, keep], D = 1)$est.theta)
  expect_equal(rst$with_purify$score, sc)
})

test_that("crdif() reports an incomplete purification when max.iter is reached", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  out <- crt_catch(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1,
                         purify = TRUE, max.iter = 1, alpha = 0.3, verbose = FALSE))
  rst <- out$value
  expect_true(any(grepl("maximum number of iterations", out$warnings)))
  expect_false(rst$with_purify$complete)
  expect_equal(rst$with_purify$n.iter, 1)
})

test_that("crdif() applies min.resp when the scores are estimated", {
  sim <- crt_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- n_resp < 4
  expect_true(any(drop_row))
  resp <- sim$resp
  resp[drop_row, ] <- NA
  out <- crt_catch(crdif(sim$x, sim$resp, group = sim$group, focal.name = 1, min.resp = 4))
  expect_true(any(grepl("fewer than 4 responses", out$warnings)))
  r2 <- suppressWarnings(crdif(sim$x, resp, group = sim$group, focal.name = 1))
  expect_identical(out$value$no_purify$dif_stat, r2$no_purify$dif_stat)
})

test_that("crdif() applies min.resp also when the scores are supplied", {
  sim <- crt_sim()
  n_resp <- rowSums(!is.na(sim$resp))
  drop_row <- n_resp < 4
  resp <- sim$resp
  resp[drop_row, ] <- NA
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  score_na <- score
  score_na[drop_row] <- NA
  out <- crt_catch(crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1, min.resp = 4))
  expect_true(any(grepl("fewer than 4 responses", out$warnings)))
  r4 <- crdif(sim$x, resp, score = score_na, group = sim$group, focal.name = 1)
  expect_identical(out$value$no_purify$dif_stat, r4$no_purify$dif_stat)
})

test_that("crdif() purification stops with a warning when every item is flagged", {
  sim <- crt_sim()
  out <- crt_catch(crdif(sim$x[4:5, ], sim$resp[, 4:5], group = sim$group, focal.name = 1,
                         purify = TRUE, alpha = 0.999, verbose = FALSE))
  expect_true(any(grepl("All items were flagged", out$warnings)))
  expect_false(out$value$with_purify$complete)
  expect_equal(out$value$with_purify$dif_item, 1:2)
  expect_equal(out$value$with_purify$n.iter, 1)
  expect_false(anyNA(out$value$with_purify$dif_stat$crdifrs))
  rst0 <- suppressWarnings(crdif(sim$x[4:5, ], sim$resp[, 4:5], group = sim$group, focal.name = 1,
                                 alpha = 0.999))
  expect_identical(out$value$no_purify, rst0$no_purify)

  # purification works when one item is left
  rst <- suppressWarnings(crdif(sim$x[3:5, ], sim$resp[, 3:5], group = sim$group, focal.name = 1,
                                purify = TRUE, alpha = 0.3, verbose = FALSE))
  expect_false(anyNA(rst$with_purify$dif_stat$crdifrs))
})

test_that("print() of crdif results leaves the significance symbols of skipped items blank", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  rst <- crdif(sim$x, sim$resp, score = score, group = sim$group, focal.name = 1, item.skip = 2)
  out <- utils::capture.output(print(rst))
  expect_false(any(grepl("?", out, fixed = TRUE)))
})

test_that("crdif() returns the number of iterations as an integer when no item is flagged", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  none <- crdif(sim$x[1:3, ], sim$resp[, 1:3], score = score, group = sim$group, focal.name = 1,
                alpha = 1e-4, purify = TRUE, verbose = FALSE)
  expect_null(none$with_purify$dif_item)
  expect_identical(none$with_purify$n.iter, 0L)
})

test_that("crdif() stops with clear errors for invalid inputs", {
  sim <- crt_sim()
  score <- suppressWarnings(est_score(sim$x, sim$resp, D = 1)$est.theta)
  run <- function(...) {
    args <- utils::modifyList(
      list(x = sim$x, data = sim$resp, score = score, group = sim$group, focal.name = 1), list(...)
    )
    do.call(crdif, args)
  }
  resp_bad <- sim$resp
  resp_bad[1, 4] <- 4
  expect_error(run(data = resp_bad), "score categories")
  resp_bad <- sim$resp
  resp_bad[is.na(resp_bad)] <- -9
  expect_error(run(data = resp_bad), "score categories")
  expect_error(run(group = sim$group[-1]), "'group'")
  expect_error(run(focal.name = 2), "focal.name")
  expect_error(run(group = replace(sim$group, 1:10, 2)), "exactly two groups")
  expect_error(run(item.skip = 8), "item.skip")
  expect_error(run(alpha = 0), "alpha")
  expect_error(run(score = score[-1]), "'score'")
  expect_error(run(purify = TRUE, max.iter = 0), "max.iter")
  expect_error(run(min.resp = -1), "min.resp")
})
