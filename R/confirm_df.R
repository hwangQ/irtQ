# Check that the item metadata is a correctly formatted data frame with all
# components required for IRT analyses, and standardize it
confirm_df <- function(x, g2na = FALSE) {
  # change all factor variables into character variables
  x <- data.frame(x, stringsAsFactors = FALSE)
  x <- purrr::modify_if(x, is.factor, as.character)
  x[, 3] <- toupper(x[, 3])
  modelGood <- all(x[, 3] %in% c("1PLM", "2PLM", "3PLM", "DRM", "GRM", "GPCM"))
  catsGood <- all(x[, 2] >= 1)
  if (!modelGood) {
    stop(paste0(
      "At least one model name is mis-specified in the model column.\n",
      "Available model names are 1PLM, 2PLM, 3PLM, DRM, GRM, and GPCM"
    ), call. = FALSE)
  }
  if (!catsGood) {
    stop(paste0(
      "At least, one score category is less than 2 under the cats column.\n",
      "Any score category should have the number greater than 1"
    ), call. = FALSE)
  }

  # add a par.3 (guessing parameter) column when it is missing and all items are 1PLM or 2PLM
  if (ncol(x[, -c(1, 2, 3)]) == 2) {
    if (all(x[, 3] %in% c("1PLM", "2PLM"))) {
      x <- data.frame(x, par.3 = NA)
    } else {
      stop("Add a par.3 column to the item metadata 'x' unless all items are 1PLM or 2PLM.", call. = FALSE)
    }
  }

  # remove the parameter columns with all NAs except par.1
  col.na.lg <- Rfast::colAll(is.na(x[, -c(1:4)]))
  if (!all(col.na.lg)) {
    x <- x[, c(!logical(4), !col.na.lg)]
  }

  # re-assign column names
  colnames(x) <- c("id", "cats", "model", paste0("par.", 1:(ncol(x) - 3)))

  # handle the g parameters for the 1PLM and 2PLM
  if (g2na) {
    # assign NAs to the par.3 column for the 1PLM, 2PLM items
    x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- NA_real_
  } else {
    # assign 0 values to the par.3 column for the 1PLM, 2PLM items
    x[x$model %in% c("1PLM", "2PLM"), "par.3"] <- 0
  }

  # consider DRM as 3PLM
  if ("DRM" %in% x$model) {
    x$model[x$model == "DRM"] <- "3PLM"
    memo <- "All 'DRM' items are considered as '3PLM' items during the item parameter estimation. \n"
    warning(memo, call. = FALSE)
  }

  # return the results
  x
}
