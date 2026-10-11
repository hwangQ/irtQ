# This function updates the ability estimates after an item is removed during the purification
# of the DIF and IPD functions. With the pattern scoring methods, the examinees are scored
# independently, so after the first iteration only the examinees who responded to the removed
# item are rescored
dif_rescore <- function(x, data, score, loc_resp, first, D, method, range,
                        norm.prior, nquad, weights, ncore, ...) {
  if (!first && ncore == 1 && method %in% c("ML", "MLF", "WL", "MAP", "EAP")) {
    # rescore the examinees who responded to the removed item
    if (length(loc_resp) > 0L) {
      score[loc_resp] <-
        est_score(
          x = x, data = data[loc_resp, , drop = FALSE], D = D, method = method,
          range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
          ncore = ncore, ...)$est.theta
    }
  } else {
    # rescore all examinees
    score <-
      est_score(
        x = x, data = data, D = D, method = method,
        range = range, norm.prior = norm.prior, nquad = nquad, weights = weights,
        ncore = ncore, ...)$est.theta
  }

  # return the updated ability estimates
  score
}
