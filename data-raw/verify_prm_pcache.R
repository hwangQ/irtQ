# =============================================================================
# Verification script for B4: nlminb P(theta) caching — PRM phase
#
# Step 2a (Numerical Equivalence):
#   For each model branch (GRM !fix.a, GPCM !fix.a, GPCM fix.a=PCM):
#   - loglike_prm, grad_item_prm, hess_item_prm with prob_cache = NULL
#     must equal those with prob_cache explicitly supplied (identical()).
#
# Step 2b (Performance):
#   microbenchmark the cached vs legacy path at three theta-grid sizes.
# =============================================================================

if (!requireNamespace("microbenchmark", quietly = TRUE)) install.packages("microbenchmark")
library(microbenchmark)

devtools::load_all(".")

# --------------------------------------------------------------------------
# Step 2a — bit-exact equivalence verification
# --------------------------------------------------------------------------

verify_prm_branch <- function(label, pr.mod, fix.a = FALSE, a.val = 1,
                              n_cats = 4, n_theta = 41, n_examinees = 500) {
  set.seed(42)
  theta <- seq(-4, 4, length.out = n_theta)

  # synthesize item parameters
  if (!fix.a) {
    a_true <- 1.4
    d_true <- sort(rnorm(n_cats - 1, 0, 1))
    item_par <- c(a_true, d_true)   # (a, d1, ..., d_{K-1})
  } else {
    a_true <- a.val
    d_true <- sort(rnorm(n_cats - 1, 0, 1))
    item_par <- d_true              # only d parameters estimated
  }

  # true category probabilities under the model
  ps_true <- prm(theta, a = a_true, d = d_true, D = 1, pr.model = pr.mod)

  # observed category frequency matrix (n_theta rows, n_cats cols)
  r_i <- ps_true * n_examinees

  nstd <- n_theta

  # --- compute prob_cache manually to test the cache code path ---
  if (pr.mod == "GRM") {
    m <- length(d_true)
    allPst <- drm(theta, a = rep(a_true, m), b = d_true, g = 0, D = 1)
    P_cache <- cbind(1, allPst) - cbind(allPst, 0)
    P_cache[P_cache > 9999999999e-10] <- 9999999999e-10
    P_cache[P_cache < 1e-10] <- 1e-10
    prob_cache <- list(allPst = allPst, P = P_cache)
  } else {
    # GPCM / PCM: same computation as make_prm_optim_fns get_prob
    d_full <- c(0, d_true)
    Da <- 1 * a_true
    theta_d <- Rfast::Outer(x = theta, y = d_full, oper = "-")
    z <- Da * theta_d
    cumsum_z <- t(Rfast::colCumSums(z))
    if (any(cumsum_z > 700)) cumsum_z <- (cumsum_z / max(cumsum_z)) * 700
    if (any(cumsum_z < -700)) cumsum_z <- -(cumsum_z / min(cumsum_z)) * 700
    numer <- exp(cumsum_z)
    denom <- Rfast::rowsums(numer)
    P_cache <- numer / denom
    prob_cache <- list(theta_d = theta_d, numer = numer, denom = denom, P = P_cache)
  }

  # loglike
  ll_null   <- loglike_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                           D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val)
  ll_cached <- loglike_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                           D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val,
                           prob_cache = prob_cache)

  # gradient
  gr_null   <- grad_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                             D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val)
  gr_cached <- grad_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                             D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val,
                             prob_cache = prob_cache)

  # hessian (adjust=FALSE for determinism)
  hs_null   <- hess_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                             D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val, adjust = FALSE)
  hs_cached <- hess_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                             D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val, adjust = FALSE,
                             prob_cache = prob_cache)

  ok_ll <- identical(ll_null, ll_cached)
  ok_gr <- identical(gr_null, gr_cached)
  ok_hs <- identical(hs_null, hs_cached)

  cat(sprintf("  [%s] ll: %s | grad: %s | hess: %s\n",
              label, ok_ll, ok_gr, ok_hs))

  invisible(all(ok_ll, ok_gr, ok_hs))
}

cat("== Step 2a: numerical equivalence (identical() bit-exact) ==\n")
res <- c(
  verify_prm_branch("GRM  4-cat !fix.a",  pr.mod = "GRM",  fix.a = FALSE, n_cats = 4),
  verify_prm_branch("GRM  5-cat !fix.a",  pr.mod = "GRM",  fix.a = FALSE, n_cats = 5),
  verify_prm_branch("GPCM 4-cat !fix.a",  pr.mod = "GPCM", fix.a = FALSE, n_cats = 4),
  verify_prm_branch("GPCM 5-cat !fix.a",  pr.mod = "GPCM", fix.a = FALSE, n_cats = 5),
  verify_prm_branch("PCM  4-cat  fix.a",  pr.mod = "GPCM", fix.a = TRUE,  a.val = 1, n_cats = 4),
  verify_prm_branch("PCM  5-cat  fix.a",  pr.mod = "GPCM", fix.a = TRUE,  a.val = 1, n_cats = 5)
)
stopifnot(all(res))
cat("All PRM branches: identical() PASSED\n\n")

# --------------------------------------------------------------------------
# Step 2b — performance: single nlminb-equivalent evaluation (obj+grad+hess)
# --------------------------------------------------------------------------

bench_prm_iter <- function(pr.mod, fix.a = FALSE, a.val = 1,
                           n_cats, n_theta, n_examinees, label) {
  set.seed(7)
  theta <- seq(-4, 4, length.out = n_theta)
  if (!fix.a) {
    a_true <- 1.4; d_true <- sort(rnorm(n_cats - 1, 0, 1))
    item_par <- c(a_true, d_true)
  } else {
    a_true <- a.val; d_true <- sort(rnorm(n_cats - 1, 0, 1))
    item_par <- d_true
  }
  ps_true <- prm(theta, a = a_true, d = d_true, D = 1, pr.model = pr.mod)
  r_i <- ps_true * n_examinees
  nstd <- n_theta

  # legacy path: each function recomputes probabilities independently
  legacy <- function() {
    o  <- loglike_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                      D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val)
    gr <- grad_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                        D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val)
    hs <- hess_item_prm(item_par, r_i = r_i, theta = theta, pr.mod = pr.mod,
                        D = 1, nstd = nstd, fix.a = fix.a, a.val = a.val, adjust = FALSE)
    list(o, gr, hs)
  }

  # cached path: single probability computation shared
  fns <- make_prm_optim_fns(
    r_i = r_i, theta = theta, pr.mod = pr.mod, D = 1, nstd = nstd,
    fix.a = fix.a, a.val = a.val,
    aprior = list(dist = "lnorm", params = c(1, 0.5)),
    bprior = list(dist = "norm", params = c(0.0, 1.0)),
    use.aprior = FALSE, use.bprior = FALSE
  )
  cached <- function() {
    o  <- fns$objective(item_par)
    gr <- fns$gradient(item_par)
    hs <- fns$hessian(item_par)
    list(o, gr, hs)
  }

  # sanity-check equivalence
  L <- legacy(); C <- cached()
  stopifnot(identical(L[[1]], C[[1]]))
  stopifnot(identical(L[[2]], C[[2]]))

  cat(sprintf("\n[%s] pr.mod=%s fix.a=%s n_cats=%d n_theta=%d n_examinees=%d\n",
              label, pr.mod, fix.a, n_cats, n_theta, n_examinees))
  print(microbenchmark(legacy = legacy(), cached = cached(),
                       times = 200L, unit = "ms"))
}

cat("== Step 2b: microbenchmark (objective + gradient + hessian per nlminb iter) ==\n")
bench_prm_iter("GRM",  FALSE, 1, n_cats=4, n_theta=41, n_examinees=500,  label="GRM  small")
bench_prm_iter("GRM",  FALSE, 1, n_cats=4, n_theta=41, n_examinees=5000, label="GRM  large")
bench_prm_iter("GPCM", FALSE, 1, n_cats=4, n_theta=41, n_examinees=500,  label="GPCM small")
bench_prm_iter("GPCM", FALSE, 1, n_cats=4, n_theta=41, n_examinees=5000, label="GPCM large")
bench_prm_iter("GPCM", TRUE,  1, n_cats=4, n_theta=41, n_examinees=500,  label="PCM  small")

cat("\nAll PRM verifications complete.\n")
