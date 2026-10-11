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
    df1 <- if (K == 2) 1 else K
    out <- rbind(out, data.frame(
      crdifr = qr$chi, crdifs = qs$chi, crdifrs = qrs$chi,
      rank.r = qr$rank, rank.s = qs$rank, rank.rs = qrs$rank,
      p.crdifr = pchisq(qr$chi, df1, lower.tail = FALSE),
      p.crdifs = pchisq(qs$chi, df1, lower.tail = FALSE),
      p.crdifrs = pchisq(qrs$chi, 2 * df1, lower.tail = FALSE),
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

test_that("crdif() purification stops when every item is flagged", {
  sim <- crt_sim()
  expect_error(
    suppressWarnings(crdif(sim$x[4:5, ], sim$resp[, 4:5], group = sim$group, focal.name = 1,
                           purify = TRUE, alpha = 0.999, verbose = FALSE)),
    "All items were flagged"
  )

  # purification works when one item is left
  rst <- suppressWarnings(crdif(sim$x[3:5, ], sim$resp[, 3:5], group = sim$group, focal.name = 1,
                                purify = TRUE, alpha = 0.3, verbose = FALSE))
  expect_false(anyNA(rst$with_purify$dif_stat$crdifrs))
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
