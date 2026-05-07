# This function computes the posterior distribution of abilities for examinees
# given their item response data, item parameters, and population distribution
#' @importFrom Rfast eachrow
posterior <- function(likehd, weights, idx.std = NULL) {
  if (is.null(idx.std)) {
    # single-group: GEMV for the denominator, eachrow for column scaling
    w <- weights[, 2]
    denom <- as.vector(likehd %*% w)
    Rfast::eachrow(likehd, w, oper = "*") / denom
  } else {
    # multi-group: same approach within each group
    ntheta <- nrow(weights[[1]])
    out <- array(0, c(nrow(likehd), ntheta))
    for (g in seq_along(weights)) {
      idx <- idx.std[[g]]
      lh_g <- likehd[idx, , drop = FALSE]
      w_g <- weights[[g]][, 2]
      denom_g <- as.vector(lh_g %*% w_g)
      out[idx, ] <- Rfast::eachrow(lh_g, w_g, oper = "*") / denom_g
    }
    out
  }
}
