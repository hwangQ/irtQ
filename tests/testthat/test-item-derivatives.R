# Analytic gradients and Hessians of the item log-likelihood (log posterior)
# compared with central finite differences of a closed-form log-likelihood

fd_grad_ref <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    e <- rep(0, length(x))
    e[j] <- h
    (f(x + e) - f(x - e)) / (2 * h)
  }, numeric(1))
}

fd_jac_ref <- function(f, x, h = 1e-5) {
  J <- sapply(seq_along(x), function(j) {
    e <- rep(0, length(x))
    e[j] <- h
    (f(x + e) - f(x - e)) / (2 * h)
  })
  J <- matrix(J, nrow = length(x))
  (J + t(J)) / 2
}

drm_ref <- function(theta, a, b, g, D) g + (1 - g) / (1 + exp(-D * a * (theta - b)))

grm_ref <- function(theta, a, d, D) {
  ps <- matrix(sapply(d, function(dk) 1 / (1 + exp(-D * a * (theta - dk)))), nrow = length(theta))
  cbind(1, ps) - cbind(ps, 0)
}

gpcm_ref <- function(theta, a, d, D) {
  lz <- t(apply(outer(theta, c(0, d), "-") * D * a, 1, cumsum))
  e <- exp(lz - apply(lz, 1, max))
  e / rowSums(e)
}

rel_err <- function(x, ref) max(abs(x - ref)) / max(1, max(abs(ref)))

theta_q <- seq(-4, 4, length.out = 21)
ap_ref <- list(dist = "lnorm", params = c(0.2, 0.5))
bp_ref <- list(dist = "norm", params = c(0.3, 1.5))
gp_ref <- list(dist = "beta", params = c(5, 16))

test_that("dichotomous item gradients and Hessians match finite differences", {
  set.seed(101)
  nq <- length(theta_q)
  cases <- list(
    list(mod = "1PLM", fix.a = FALSE, fix.g = FALSE, par = c(1.2, -0.5, 0.3, 1.1), n1 = 3),
    list(mod = "1PLM", fix.a = TRUE, fix.g = FALSE, par = 0.4, a.val = 1.3),
    list(mod = "2PLM", fix.a = TRUE, fix.g = FALSE, par = c(1.4, 0.2)),
    list(mod = "3PLM", fix.a = TRUE, fix.g = FALSE, par = c(1.6, 0.3, 0.17)),
    list(mod = "3PLM", fix.a = TRUE, fix.g = TRUE, par = c(1.6, 0.3), g.val = 0.15)
  )
  for (D in c(1, 1.702)) {
    for (pri in c(FALSE, TRUE)) {
      for (cs in cases) {
        n1 <- if (is.null(cs$n1)) 1 else cs$n1
        f_i <- matrix(round(runif(nq * n1, 20, 200)), nq, n1)
        r_i <- matrix(rbinom(nq * n1, f_i, 0.5), nq, n1)
        s_i <- f_i - r_i
        if (n1 == 1) {
          f_i <- c(f_i)
          r_i <- c(r_i)
          s_i <- c(s_i)
        }
        a.val <- if (is.null(cs$a.val)) 1 else cs$a.val
        g.val <- if (is.null(cs$g.val)) 0.2 else cs$g.val
        common <- list(
          r_i = r_i, theta = theta_q, mod = cs$mod, D = D, fix.a = cs$fix.a,
          fix.g = cs$fix.g, a.val = a.val, g.val = g.val, n.1PLM = cs$n1,
          aprior = ap_ref, bprior = bp_ref, gprior = gp_ref,
          use.aprior = pri, use.bprior = pri, use.gprior = pri
        )
        gr <- function(p) do.call(irtQ:::grad_item_drm, c(list(item_par = p, f_i = f_i), common))
        he <- function(p) do.call(irtQ:::hess_item_drm_inner, c(list(item_par = p, f_i = f_i), common))
        nll <- function(p) do.call(irtQ:::loglike_drm, c(list(item_par = p, s_i = s_i), common))

        # closed-form negative log posterior
        ref <- function(p) {
          if (cs$mod == "1PLM" && !cs$fix.a) {
            P <- sapply(p[-1], function(b) drm_ref(theta_q, p[1], b, 0, D))
            ll <- sum(r_i * log(P) + s_i * log(1 - P))
            if (pri) {
              ll <- ll + dlnorm(D * p[1], 0.2, 0.5, log = TRUE) + sum(dnorm(p[-1], 0.3, 1.5, log = TRUE))
            }
            return(-ll)
          }
          a <- switch(cs$mod, "1PLM" = a.val, p[1])
          b <- switch(cs$mod, "1PLM" = p[1], p[2])
          g <- if (cs$mod != "3PLM") 0 else if (cs$fix.g) g.val else p[3]
          P <- drm_ref(theta_q, a, b, g, D)
          ll <- sum(r_i * log(P) + s_i * log(1 - P))
          if (pri) {
            if (cs$mod != "1PLM") ll <- ll + dlnorm(D * a, 0.2, 0.5, log = TRUE)
            ll <- ll + dnorm(b, 0.3, 1.5, log = TRUE)
            if (cs$mod == "3PLM" && !cs$fix.g) ll <- ll + dbeta(g, 5, 16, log = TRUE)
          }
          -ll
        }
        p0 <- cs$par
        expect_equal(nll(p0), ref(p0), tolerance = 1e-10)
        expect_lt(rel_err(gr(p0), fd_grad_ref(ref, p0)), 1e-7)
        expect_lt(rel_err(he(p0), fd_jac_ref(gr, p0)), 1e-7)
      }
    }
  }
})

test_that("polytomous item gradients and Hessians match finite differences", {
  set.seed(102)
  nq <- length(theta_q)
  cases <- list(
    list(mod = "GRM", K = 3, par = c(1.2, -0.7, 0.8)),
    list(mod = "GRM", K = 5, par = c(1.5, -1.5, -0.4, 0.5, 1.7)),
    list(mod = "GPCM", K = 3, par = c(1.2, -0.7, 0.8)),
    list(mod = "GPCM", K = 5, par = c(1.3, -1.2, 0.3, -0.2, 1.1)),
    list(mod = "GPCM", K = 4, par = c(-0.5, 0.2, 0.9), fix.a = TRUE, a.val = 1.4),
    list(mod = "GPCM", K = 2, par = 0.3, fix.a = TRUE, a.val = 1)
  )
  for (D in c(1, 1.702)) {
    for (pri in c(FALSE, TRUE)) {
      for (cs in cases) {
        K <- cs$K
        fa <- isTRUE(cs$fix.a)
        a.val <- if (fa) cs$a.val else 1
        r_i <- matrix(round(runif(nq * K, 5, 80)), nq, K)
        common <- list(
          r_i = r_i, theta = theta_q, pr.mod = cs$mod, D = D, fix.a = fa, a.val = a.val,
          aprior = ap_ref, bprior = bp_ref, use.aprior = pri, use.bprior = pri
        )
        nll <- function(p) do.call(irtQ:::loglike_prm, c(list(item_par = p), common))
        gr <- function(p) do.call(irtQ:::grad_item_prm, c(list(item_par = p), common))
        he <- function(p) do.call(irtQ:::hess_item_prm_inner, c(list(item_par = p, nstd = nq), common))
        ref <- function(p) {
          a <- if (fa) a.val else p[1]
          d <- if (fa) p else p[-1]
          P <- if (cs$mod == "GRM") grm_ref(theta_q, a, d, D) else gpcm_ref(theta_q, a, d, D)
          ll <- sum(r_i * log(P))
          if (pri) {
            if (!fa) ll <- ll + dlnorm(D * a, 0.2, 0.5, log = TRUE)
            ll <- ll + sum(dnorm(d, 0.3, 1.5, log = TRUE))
          }
          -ll
        }
        p0 <- cs$par
        expect_equal(nll(p0), ref(p0), tolerance = 1e-10)
        expect_lt(rel_err(gr(p0), fd_grad_ref(ref, p0)), 1e-7)
        expect_lt(rel_err(he(p0), fd_jac_ref(gr, p0)), 1e-7)
      }
    }
  }
})

test_that("GPCM derivatives stay finite and consistent for large slopes", {
  th <- seq(-6, 6, length.out = 49)
  r_i <- matrix(1, 49, 4)
  ref <- function(p) -sum(r_i * log(pmax(gpcm_ref(th, p[1], p[-1], 1.702), .Machine$double.xmin)))
  for (a in c(5, 8, 20, 100)) {
    p0 <- c(a, -1, 0, 1)
    fns <- irtQ:::make_prm_optim_fns(
      r_i = r_i, theta = th, pr.mod = "GPCM", D = 1.702, nstd = 49, fix.a = FALSE, a.val = 1,
      aprior = ap_ref, bprior = bp_ref, use.aprior = FALSE, use.bprior = FALSE
    )
    # the cached objective equals the objective computed from prm()
    expect_equal(fns$objective(p0), irtQ:::loglike_prm(p0, r_i, th, "GPCM", D = 1.702), tolerance = 1e-10)
    expect_equal(fns$objective(p0), ref(p0), tolerance = 1e-10)
    g <- fns$gradient(p0)
    expect_true(all(is.finite(g)))
    expect_lt(rel_err(g, fd_grad_ref(ref, p0, 1e-7)), 1e-6)

    # the Hessian is finite while no category probability underflows
    if (a <= 8) {
      h <- irtQ:::hess_item_prm_inner(p0, r_i, th, "GPCM", D = 1.702, nstd = 49)
      expect_true(all(is.finite(h)))
      expect_lt(rel_err(h, fd_jac_ref(fns$gradient, p0, 1e-6)), 1e-6)
    }
  }
})
