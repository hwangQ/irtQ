# ---- Fixtures ------------------------------------------------------------------

prm_file <- system.file("extdata", "flexmirt_sample-prm.txt", package = "irtQ")
x_full   <- bring.flexmirt(file = prm_file, "par")$Group1$full_df  # 55 items

## 1. 3PLM-only: 10 dichotomous items (from flexMIRT sample)
x_drm <- x_full[1:10, ]
set.seed(99)
theta_drm  <- rnorm(200)
resp_drm   <- simdat(x = x_drm, theta = theta_drm, D = 1)

## 2. GRM-only: 5 items, 4 categories (scored 0-3)
x_grm <- shape_df(
  par.prm = list(
    a = c(0.9, 1.1, 1.2, 1.0, 0.8),
    d = list(c(-1.2, -0.2, 0.8), c(-1.0, 0.0, 1.0),
             c(-0.8,  0.2, 1.2), c(-1.5,-0.3, 0.7),
             c(-0.9,  0.1, 0.9))
  ),
  cats = rep(4L, 5), model = "GRM"
)
set.seed(11)
theta_grm  <- rnorm(300)
resp_grm   <- simdat(x = x_grm, theta = theta_grm, D = 1)

## 3. GPCM-only: 5 items, 4 categories
x_gpcm <- shape_df(
  par.prm = list(
    a = c(1.0, 1.2, 0.9, 1.1, 1.0),
    d = list(c(-1.0, 0.0, 0.8), c(-0.8, 0.2, 1.0),
             c(-1.2,-0.2, 0.9), c(-0.5, 0.5, 1.2),
             c(-1.0,-0.1, 0.7))
  ),
  cats = rep(4L, 5), model = "GPCM"
)
set.seed(22)
theta_gpcm <- rnorm(300)
resp_gpcm  <- simdat(x = x_gpcm, theta = theta_gpcm, D = 1)

## 4. Mixed: 5 x 3PLM + 2 x GRM (4 cats)  -> max sum score = 5 + 6 = 11
x_mix_grm <- shape_df(
  par.drm = list(a = c(1.0,1.2,0.9,1.1,1.0),
                 b = c(-1.0,-0.5,0.0,0.5,1.0),
                 g = c(0.2, 0.15,0.1,0.2,0.15)),
  par.prm = list(a = c(1.0,1.2),
                 d = list(c(-1.0,0.0,1.0), c(-0.5,0.5,1.5))),
  cats  = c(2L,2L,2L,2L,2L,4L,4L),
  model = c(rep("3PLM",5), "GRM","GRM")
)
set.seed(33)
theta_mix_grm  <- rnorm(300)
resp_mix_grm   <- simdat(x = x_mix_grm, theta = theta_mix_grm, D = 1)

## 5. Mixed: 5 x 3PLM + 2 x GPCM (4 cats) -> max sum score = 5 + 6 = 11
x_mix_gpcm <- shape_df(
  par.drm = list(a = c(1.0,1.2,0.9,1.1,1.0),
                 b = c(-1.0,-0.5,0.0,0.5,1.0),
                 g = c(0.2, 0.15,0.1,0.2,0.15)),
  par.prm = list(a = c(1.0,1.2),
                 d = list(c(-1.0,0.0,1.0), c(-0.5,0.5,1.5))),
  cats  = c(2L,2L,2L,2L,2L,4L,4L),
  model = c(rep("3PLM",5), "GPCM","GPCM")
)
set.seed(44)
theta_mix_gpcm  <- rnorm(300)
resp_mix_gpcm   <- simdat(x = x_mix_gpcm, theta = theta_mix_gpcm, D = 1)


# ---- Helpers -------------------------------------------------------------------

check_pointwise <- function(result, n = 200) {
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n)
  expect_true("est.theta" %in% colnames(result))
  expect_true("se.theta"  %in% colnames(result))
  expect_true(all(!is.nan(result$se.theta) | is.na(result$se.theta)))
}

check_sumtable <- function(result, max_ss) {
  expect_type(result, "list")
  expect_true("est.par"     %in% names(result))
  expect_true("score.table" %in% names(result))
  expect_true(all(c("sum.score","est.theta","se.theta") %in% colnames(result$est.par)))
  expect_equal(nrow(result$score.table), max_ss + 1L)
  expect_true(all(diff(result$score.table$est.theta) >= 0))
}


# ==============================================================================
# 1. 3PLM ONLY (10 dichotomous items, max sum = 10)
# ==============================================================================

base_drm <- list(x = x_drm, data = resp_drm, D = 1)

test_that("3PLM | ML: returns est.theta/se.theta data frame, correlates with truth", {
  res <- do.call(est_score, c(base_drm, list(method = "ML", range = c(-6, 6))))
  check_pointwise(res)
  expect_gt(cor(res$est.theta, theta_drm), 0.7)
})

test_that("3PLM | MLF: finite estimates within fence range", {
  res <- do.call(est_score, c(base_drm, list(method = "MLF", range = c(-5, 5))))
  check_pointwise(res)
  expect_true(all(is.finite(res$est.theta)))
  expect_true(all(res$est.theta >= -5 & res$est.theta <= 5))
})

test_that("3PLM | WL: returns est.theta/se.theta data frame", {
  res <- do.call(est_score, c(base_drm, list(method = "WL", range = c(-6, 6))))
  check_pointwise(res)
})

test_that("3PLM | MAP: returns est.theta/se.theta data frame", {
  res <- do.call(est_score, c(base_drm, list(method = "MAP",
                                              norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res)
})

test_that("3PLM | EAP: estimates within quadrature range", {
  res <- do.call(est_score, c(base_drm, list(method = "EAP",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res)
  expect_true(all(res$est.theta >= -4 & res$est.theta <= 4))
})

test_that("3PLM | EAP.SUM: list structure, score.table has 11 rows (0-10)", {
  res <- do.call(est_score, c(base_drm, list(method = "EAP.SUM",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 10L)
  expect_equal(nrow(res$est.par), 200L)
})

test_that("3PLM | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_drm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 10L)
})

test_that("3PLM | ML handles partial NA responses", {
  resp_na <- resp_drm; resp_na[1,1] <- NA; resp_na[2,2:3] <- NA
  res <- est_score(x=x_drm, data=resp_na, D=1, method="ML", missing=NA, range=c(-6,6))
  check_pointwise(res)
})

test_that("3PLM | EAP handles partial NA responses", {
  resp_na <- resp_drm; resp_na[5, 1:3] <- NA
  res <- est_score(x=x_drm, data=resp_na, D=1, method="EAP",
                   missing=NA, norm.prior=c(0,1))
  check_pointwise(res)
})

test_that("3PLM | ML/WL/MAP/EAP rank-correlate > 0.95", {
  ml  <- est_score(x=x_drm, data=resp_drm, D=1, method="ML",  range=c(-6,6))$est.theta
  wl  <- est_score(x=x_drm, data=resp_drm, D=1, method="WL",  range=c(-6,6))$est.theta
  map <- est_score(x=x_drm, data=resp_drm, D=1, method="MAP", range=c(-6,6),
                   norm.prior=c(0,1))$est.theta
  eap <- est_score(x=x_drm, data=resp_drm, D=1, method="EAP",
                   norm.prior=c(0,1))$est.theta
  expect_gt(cor(ml, wl,  method="spearman"), 0.95)
  expect_gt(cor(ml, map, method="spearman"), 0.95)
  expect_gt(cor(ml, eap, method="spearman"), 0.95)
})


# ==============================================================================
# 2. GRM ONLY (5 items, 4 cats -> max sum = 15)
# ==============================================================================

base_grm <- list(x = x_grm, data = resp_grm, D = 1)

test_that("GRM | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_grm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "MAP",
                                              norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GRM | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_grm, list(method = "EAP",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("GRM | EAP.SUM: list structure, score.table has 16 rows (0-15)", {
  res <- do.call(est_score, c(base_grm, list(method = "EAP.SUM",
                                              norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 15L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("GRM | INV.TCC: list structure, score.table monotone (16 rows)", {
  res <- do.call(est_score, c(base_grm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 15L)
})

test_that("GRM | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_grm, list(method = "ML", range = c(-6,6))))
  # 5 polytomous items: lower information than 10 dichotomous -> threshold 0.65
  expect_gt(cor(res$est.theta, theta_grm), 0.65)
})


# ==============================================================================
# 3. GPCM ONLY (5 items, 4 cats -> max sum = 15)
# ==============================================================================

base_gpcm <- list(x = x_gpcm, data = resp_gpcm, D = 1)

test_that("GPCM | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "MAP",
                                               norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_gpcm, list(method = "EAP",
                                               norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("GPCM | EAP.SUM: list structure, score.table has 16 rows (0-15)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "EAP.SUM",
                                               norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 15L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("GPCM | INV.TCC: list structure, score.table monotone (16 rows)", {
  res <- do.call(est_score, c(base_gpcm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 15L)
})

test_that("GPCM | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_gpcm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_gpcm), 0.7)
})


# ==============================================================================
# 4. MIXED: 5 x 3PLM + 2 x GRM (max sum = 5*1 + 2*3 = 11)
# ==============================================================================

base_mix_grm <- list(x = x_mix_grm, data = resp_mix_grm, D = 1)

test_that("Mixed(3PLM+GRM) | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "MAP",
                                                  norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "EAP",
                                                  norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GRM) | EAP.SUM: list structure, score.table has 12 rows (0-11)", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "EAP.SUM",
                                                  norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 11L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("Mixed(3PLM+GRM) | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 11L)
})

test_that("Mixed(3PLM+GRM) | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_mix_grm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_mix_grm), 0.7)
})


# ==============================================================================
# 5. MIXED: 5 x 3PLM + 2 x GPCM (max sum = 5*1 + 2*3 = 11)
# ==============================================================================

base_mix_gpcm <- list(x = x_mix_gpcm, data = resp_mix_gpcm, D = 1)

test_that("Mixed(3PLM+GPCM) | ML: returns correct data frame structure (n=300)", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "ML", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | WL: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "WL", range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | MAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "MAP",
                                                   norm.prior = c(0,1), range = c(-6,6))))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | EAP: returns correct data frame structure", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "EAP",
                                                   norm.prior = c(0,1), nquad = 41L)))
  check_pointwise(res, n = 300L)
})

test_that("Mixed(3PLM+GPCM) | EAP.SUM: list structure, score.table has 12 rows (0-11)", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "EAP.SUM",
                                                   norm.prior = c(0,1), nquad = 41L)))
  check_sumtable(res, max_ss = 11L)
  expect_equal(nrow(res$est.par), 300L)
})

test_that("Mixed(3PLM+GPCM) | INV.TCC: list structure, score.table monotone", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "INV.TCC", range.tcc = c(-7,7))))
  check_sumtable(res, max_ss = 11L)
})

test_that("Mixed(3PLM+GPCM) | ML estimates correlate with true theta", {
  res <- do.call(est_score, c(base_mix_gpcm, list(method = "ML", range = c(-6,6))))
  expect_gt(cor(res$est.theta, theta_mix_gpcm), 0.7)
})


# ==============================================================================
# 6. INPUT CHECKS AND PARALLEL SCORING
# ==============================================================================

test_that("est_score() stops when data and item metadata have different numbers of items", {
  for (m in c("ML", "EAP", "EAP.SUM", "INV.TCC")) {
    expect_error(est_score(x_drm, resp_drm[1:3, 1:9], D = 1, method = m), "number of columns")
    expect_error(est_score(x_drm, cbind(resp_drm[1:3, ], 1), D = 1, method = m), "number of columns")
  }

  # a single response vector is checked in the same way
  expect_error(est_score(x_drm, resp_drm[1, 1:5], D = 1), "number of columns")
})

test_that("est_score() stops for responses that are not integer scores within the categories", {
  r <- resp_drm[1:3, ]

  # a score above the highest category of a dichotomous item
  r_high <- r
  r_high[1, 2] <- 2
  expect_error(est_score(x_drm, r_high, D = 1), "Responses must be integers")
  expect_error(est_score(x_drm, r_high, D = 1, method = "EAP"), "Responses must be integers")
  expect_error(est_score(x_drm, r_high, D = 1, method = "EAP.SUM"), "Responses must be integers")

  # a non-integer score
  r_frac <- r
  r_frac[2, 3] <- 0.7
  expect_error(est_score(x_drm, r_frac, D = 1), "Responses must be integers")

  # a negative code that is not declared as missing
  r_neg <- r
  r_neg[3, 1] <- -9
  expect_error(est_score(x_drm, r_neg, D = 1), "Responses must be integers")

  # the same code is accepted when it is declared as missing
  expect_identical(
    est_score(x_drm, r_neg, D = 1, missing = -9),
    est_score(x_drm, replace(r_neg, r_neg == -9, NA), D = 1)
  )

  # a polytomous item accepts its own categories only
  r_grm <- resp_grm[1:3, ]
  expect_no_error(est_score(x_grm, r_grm, D = 1))
  r_grm[1, 1] <- 4
  expect_error(est_score(x_grm, r_grm, D = 1), "Responses must be integers")
})

test_that("est_score() stops for an unknown method or starting value option", {
  r <- resp_drm[1:3, ]
  expect_error(est_score(x_drm, r, D = 1, method = "ml"), "'method' must be one of")
  expect_error(est_score(x_drm, r, D = 1, method = c("ML", "WL")), "'method' must be one of")
  expect_error(est_score(x_drm, r, D = 1, stval.opt = 4), "'stval.opt' must be 1, 2, or 3")
  expect_error(est_score(x_drm, r, D = 1, stval.opt = c(1, 2)), "'stval.opt' must be 1, 2, or 3")
})

test_that("est_score() for an est_irt object uses the stored data and D and rejects data or D", {
  fit <- structure(list(par.est = x_drm, data = resp_drm[1:5, ], scale.D = 1.702), class = "est_irt")
  expect_identical(
    est_score(fit, method = "ML"),
    est_score(x_drm, resp_drm[1:5, ], D = 1.702, method = "ML")
  )
  expect_error(est_score(fit, data = resp_drm[1:2, ], method = "ML"), "cannot be supplied")
  expect_error(est_score(fit, D = 1, method = "ML"), "cannot be supplied")
})

test_that("parallel scoring gives the serial result for any number of examinees per core", {
  skip_on_cran()
  n_conn <- nrow(showConnections())
  for (n in c(1, 2, 3, 5)) {
    r <- resp_drm[seq_len(n), , drop = FALSE]
    serial <- est_score(x_drm, r, D = 1, method = "ML")
    for (nc in c(2, 4)) {
      expect_warning(par <- est_score(x_drm, r, D = 1, method = "ML", ncore = nc), "ncore > 1")
      expect_equal(par, serial)
    }
  }

  # the clusters are closed after scoring
  expect_equal(nrow(showConnections()), n_conn)
})

test_that("lwrc() accepts a single theta for a test with polytomous items", {
  x_pm <- x_full[c(1, 53), ]
  one <- lwrc(x_pm, theta = 0.4)
  two <- lwrc(x_pm, theta = c(0.4, 1))
  expect_equal(dim(one), c(6L, 1L))
  expect_equal(unname(one[, 1]), unname(two[, 1]))
  expect_equal(sum(one), 1)
})

test_that("llike_score() uses default fences for MLF and keeps examinees with duplicated row names apart", {
  x_ll <- x_full[c(1:6, 51:53), ]
  set.seed(8)
  d_ll <- simdat(x_ll, rnorm(3), D = 1)
  th <- seq(-2, 2, 1)

  # the default fences are c(-5, 5)
  expect_identical(
    llike_score(x_ll, d_ll, th, method = "MLF"),
    llike_score(x_ll, d_ll, th, method = "MLF", fence.b = c(-5, 5))
  )

  # two identical response rows with the same row name give two identical columns
  d_dup <- rbind(d_ll[1, ], d_ll[1, ], d_ll[2, ])
  rownames(d_dup) <- c("a", "a", "b")
  ll_dup <- llike_score(x_ll, d_dup, th)
  ref <- llike_score(x_ll, d_ll[1:2, ], th)
  expect_equal(unname(as.matrix(ll_dup)), unname(as.matrix(ref[, c(1, 1, 2)])))
})

test_that("llike_score() returns NA for an examinee without any observed response", {
  x_ll <- x_full[c(1:6, 51:53), ]
  set.seed(8)
  d_ll <- simdat(x_ll, rnorm(3), D = 1)
  th <- seq(-2, 2, 1)
  d_na <- rbind(d_ll[1, ], rep(NA, ncol(d_ll)), d_ll[2, ])
  for (m in c("ML", "MAP")) {
    res <- llike_score(x_ll, d_na, th, method = m)
    expect_true(all(is.na(res[[2]])))
    expect_equal(res[, c(1, 3)], llike_score(x_ll, d_ll[1:2, ], th, method = m), ignore_attr = TRUE)
  }
})
