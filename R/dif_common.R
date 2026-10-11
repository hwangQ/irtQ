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


# This function stops when the group vector or the focal group is invalid for a two-group
# analysis of the DIF and IPD functions
check_group <- function(group, focal.name, nrow_data) {
  # stop when the group vector does not match the rows of the response data
  if (length(group) != nrow_data) {
    stop("The length of 'group' must equal the number of rows in 'data'.", call. = FALSE)
  }

  # stop when the focal group is not a single value found in the group vector
  if (length(focal.name) != 1L || !any(group == focal.name, na.rm = TRUE)) {
    stop("'focal.name' must be a single value found in 'group'.", call. = FALSE)
  }

  # stop when there is no examinee in the reference group
  if (!any(group != focal.name, na.rm = TRUE)) {
    stop("'group' must contain at least one examinee of the reference group.", call. = FALSE)
  }
}


# This function drops the item parameter columns that no remaining item uses, which are the
# columns beyond the larger of three and the maximum number of score categories
trim_par_cols <- function(x) {
  # find the parameter columns that are not used by any item
  n_par <- max(3, max(x$cats))
  par_col <- grep("^par\\.", names(x))
  drop_col <- par_col[seq_along(par_col) > n_par]

  # remove the unused columns
  if (length(drop_col) > 0L) {
    x <- x[, -drop_col, drop = FALSE]
  }
  x
}


# This function checks the inputs of rdif() and crdif(), recodes the missing responses, and
# returns the response matrix, the ability estimates, and the positions of the items to be skipped
dif_prepare <- function(x, data, score, group, focal.name, item.skip, D, alpha, missing,
                        purify, max.iter, min.resp, method, range, norm.prior, nquad,
                        weights, ncore, ...) {
  # transform the response data to a matrix form
  data <- as.matrix(data)

  # re-code missing values
  if (!is.na(missing)) {
    data[data == missing] <- NA
  }

  # transform the response data to a numeric matrix and check the responses
  data <- resp_to_matrix(data, x$cats, x$id)

  # check the group vector and the focal group
  check_group(group, focal.name, nrow(data))

  # stop when the group vector has more than two groups
  if (length(unique(group[!is.na(group)])) > 2L) {
    stop("'group' must contain exactly two groups (the reference and focal groups). ",
         "Use grdif() for more than two groups.", call. = FALSE)
  }

  # stop when the significance level is not a single value between 0 and 1
  if (!is.numeric(alpha) || length(alpha) != 1L || is.na(alpha) || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number between 0 and 1.", call. = FALSE)
  }

  # stop when the maximum number of iterations is not a positive whole number
  if (purify && !(is.numeric(max.iter) && length(max.iter) == 1L && is.finite(max.iter) &&
                  max.iter >= 1 && max.iter == round(max.iter))) {
    stop("'max.iter' must be a single positive whole number when purify = TRUE.", call. = FALSE)
  }

  # stop when the minimum number of responses is not a single nonnegative number
  if (!is.null(min.resp) && !(is.numeric(min.resp) && length(min.resp) == 1L &&
                              is.finite(min.resp) && min.resp >= 0)) {
    stop("'min.resp' must be NULL or a single nonnegative number.", call. = FALSE)
  }

  # check the positions of the items to be skipped
  item.skip <- check_item_skip(item.skip, nrow(x))

  # compute the score if score = NULL
  if (!is.null(score)) {
    # transform scores to a vector form
    if (is.matrix(score) | is.data.frame(score)) {
      score <- as.numeric(data.matrix(score))
    }

    # stop when the scores are not numbers
    if (!is.numeric(score)) {
      stop("'score' must be a numeric vector.", call. = FALSE)
    }

    # stop when the scores do not match the rows of the response data
    if (length(score) != nrow(data)) {
      stop("The length of 'score' must equal the number of rows in 'data'.", call. = FALSE)
    }
  } else {
    # if min.resp is not NULL, find the examinees who have the number of responses
    # less than specified value (e.g., 5). Then, replace their all responses with NA
    if (!is.null(min.resp)) {
      n_resp <- rowSums(!is.na(data))
      loc_less <- which(n_resp < min.resp & n_resp > 0)

      # warn that the examinees with too few responses are excluded
      if (length(loc_less) > 0L) {
        warning(length(loc_less), " examinee(s) with fewer than ", min.resp,
                " responses were excluded.", call. = FALSE)
      }
      data[loc_less, ] <- NA
    }
    score <- est_score(
      x = x, data = data, D = D, method = method, range = range, norm.prior = norm.prior,
      nquad = nquad, weights = weights, ncore = ncore, ...)$est.theta
  }

  # stop when no examinee has an ability estimate
  if (all(is.na(score[!is.na(group)]))) {
    stop("No examinee has an ability estimate.", call. = FALSE)
  }

  # exclude the examinees with responses but without an ability estimate, with a warning
  na_score <- is.na(score) & rowSums(!is.na(data)) > 0
  if (any(na_score)) {
    warning(sum(na_score), " examinee(s) with a missing ability estimate were excluded.",
            call. = FALSE)
    data[na_score, ] <- NA
  }

  # return the checked inputs
  list(data = data, score = score, item.skip = item.skip)
}


# This function stops when the response data or the scaling factor D is passed through '...'
# to a method for an estimation object, which takes them from the object
check_obj_dots <- function(dots) {
  if (any(c("data", "D") %in% names(dots))) {
    stop("The response data and the scaling factor D are taken from the object 'x'.",
         call. = FALSE)
  }
}


# This function computes the quadratic form of a deviation vector in the inverse of a covariance
# matrix. When the matrix is singular or nearly singular, it uses the Moore-Penrose generalized
# inverse and returns the rank of the matrix as the degrees of freedom
quad_form <- function(cov_mat, dev_vec, df) {
  # return a missing statistic when the covariance matrix has non-finite values
  if (!all(is.finite(cov_mat))) {
    return(list(stat = NA_real_, df = df, reduced = FALSE))
  }

  # compute the reciprocal condition number of the covariance matrix
  rc_cov <- tryCatch(rcond(cov_mat), error = function(e) 0)

  # use the inverse of the covariance matrix when it is well conditioned
  if (rc_cov > 1e-10) {
    stat <- as.numeric(t(dev_vec) %*% solve(cov_mat) %*% dev_vec)
    return(list(stat = stat, df = df, reduced = FALSE))
  }

  # otherwise use the generalized inverse and the rank of the covariance matrix
  eig <- eigen(cov_mat, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-10
  if (!any(keep)) {
    return(list(stat = NA_real_, df = df, reduced = FALSE))
  }
  proj <- crossprod(eig$vectors[, keep, drop = FALSE], dev_vec)
  list(stat = sum(proj^2 / eig$values[keep]), df = sum(keep), reduced = sum(keep) < df)
}
