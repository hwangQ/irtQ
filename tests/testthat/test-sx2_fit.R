# sx2_fit() computes the S-X2 item fit statistic of Orlando and Thissen (2000)

x_lsat <- shape_df(
  par.drm = list(a = c(0.8, 0.7, 0.9, 0.6, 0.7), b = c(-3.4, -1.3, -0.3, -1.8, -2.8), g = rep(0, 5)),
  item.id = paste0("V", 1:5), cats = 2, model = "2PLM"
)

test_that("sx2_fit() stops for data that do not match the items or the score categories", {
  # a column is missing
  expect_error(sx2_fit(x_lsat, data = LSAT6[, 1:4], D = 1), "number of columns")

  # a response above the highest category of an item
  resp <- LSAT6
  resp[, 1] <- resp[, 1] + 1
  expect_error(sx2_fit(x_lsat, data = resp, D = 1), "outside the score categories")

  # a non-integer response
  resp <- LSAT6
  resp[3, 2] <- 0.5
  expect_error(sx2_fit(x_lsat, data = resp, D = 1), "outside the score categories")
})

test_that("sx2_fit() reads character and factor responses like numeric responses", {
  resp_chr <- matrix(as.character(as.matrix(LSAT6)), nrow(LSAT6))
  df_chr <- as.data.frame(resp_chr, stringsAsFactors = FALSE)
  df_fac <- as.data.frame(lapply(df_chr, factor))
  ref <- sx2_fit(x_lsat, data = LSAT6, D = 1)
  for (dat in list(resp_chr, df_chr, df_fac)) {
    fit <- sx2_fit(x_lsat, data = dat, D = 1)
    expect_identical(fit$fit_stat, ref$fit_stat)
    expect_identical(unname(fit$exp_freq), unname(ref$exp_freq))
    expect_identical(unname(fit$obs_freq), unname(ref$obs_freq))
  }

  # a response that cannot be read as a number stops
  resp_chr[1, 1] <- "a"
  expect_error(sx2_fit(x_lsat, data = resp_chr, D = 1), "numeric scores")
})

test_that("sx2_fit() replaces missing responses with zeros and names the items without any response", {
  resp <- LSAT6
  resp[1:5, 2] <- NA
  resp0 <- resp
  resp0[is.na(resp0)] <- 0
  expect_warning(f1 <- sx2_fit(x_lsat, data = resp, D = 1), "Missing responses are replaced with 0")
  expect_identical(f1$fit_stat, sx2_fit(x_lsat, data = resp0, D = 1)$fit_stat)

  # an item without any response is named in the warning
  resp[, 4] <- NA
  expect_warning(sx2_fit(x_lsat, data = resp, D = 1), "item\\(s\\) V4")
})
