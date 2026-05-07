# This function computes a matrix of likelihoods for obtaining the observed item
# responses across examinees given theta values, used by the MMLE-EM algorithm.
# Returns list(L = ...) — the LL slot was removed because no caller reads it.
#' @importFrom Matrix tcrossprod
likelihood <- function(elm_item, idx.drm = NULL, idx.prm = NULL,
                       data_drm = NULL, data_prm = NULL, theta, D = 1) {
  # log-likelihood contribution from DRM items
  if (!is.null(idx.drm)) {
    ps <- drm(
      theta = theta, a = elm_item$pars[idx.drm, 1], b = elm_item$par[idx.drm, 2],
      g = elm_item$par[idx.drm, 3], D = D
    )
    qs <- 1 - ps
    log_qps <- log(cbind(qs, ps))
    llike_drm <- Matrix::tcrossprod(x = data_drm, y = log_qps)
  } else {
    llike_drm <- 0L
  }

  # log-likelihood contribution from PRM items
  if (!is.null(idx.prm)) {
    n.prm <- length(idx.prm)
    prob.prm <- vector("list", n.prm)
    for (k in 1:n.prm) {
      par.tmp <- stats::na.exclude(elm_item$par[idx.prm[k], ])
      prob.prm[[k]] <- prm(
        theta = theta, a = par.tmp[1], d = par.tmp[-1],
        D = D, pr.model = elm_item$model[idx.prm[k]]
      )
    }
    prob.prm <- do.call(what = "cbind", prob.prm)
    log_prob.prm <- log(prob.prm)
    llike_prm <- Matrix::tcrossprod(x = data_prm, y = log_prob.prm)
  } else {
    llike_prm <- 0L
  }

  # combine and exponentiate in a single pass; LL itself is no longer needed
  L <- exp(base::as.matrix(llike_drm + llike_prm))

  list(L = L)
}
